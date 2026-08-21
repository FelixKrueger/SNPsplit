//! Locating N-masked positions in a read, from its MD tag and its CIGAR.
//!
//! The MD tag records mismatches against the reference, and an N-masked genome turns every
//! known SNP into a mismatch spelled `N`. Finding where those Ns fall in the read, and what
//! genomic position each corresponds to, is the whole of the allele call.
//!
//! Two things make this harder than it looks, and both are why the Perl has two code paths:
//! insertions and soft-clips consume read bases but do not appear in the MD tag at all, so a
//! position derived from the MD tag alone points at the wrong base; and deletions appear in
//! the MD tag but do not consume read bases.

/// One CIGAR operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Op {
    pub kind: u8,
    pub len: usize,
}

/// Split a CIGAR string into operations.
///
/// Returns `None` if the lengths and operations do not pair up, which the Perl treats as
/// fatal.
pub fn parse_cigar(cigar: &str) -> Option<Vec<Op>> {
    let mut ops = Vec::new();
    let mut digits = String::new();
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else {
            if digits.is_empty() {
                return None;
            }
            ops.push(Op {
                kind: c as u8,
                len: digits.parse().ok()?,
            });
            digits.clear();
        }
    }
    digits.is_empty().then_some(ops)
}

/// Perl's `split /PATTERN/` semantics: trailing empty fields are dropped.
///
/// This is not a detail. `split /[ATCG]/, "5A"` yields one field in Perl and two in Rust, and
/// the count feeds directly into the position arithmetic below.
fn perl_split(text: &str, is_separator: impl Fn(char) -> bool) -> Vec<&str> {
    let mut parts: Vec<&str> = Vec::new();
    let mut start = 0;
    for (index, c) in text.char_indices() {
        if is_separator(c) {
            parts.push(&text[start..index]);
            start = index + c.len_utf8();
        }
    }
    parts.push(&text[start..]);
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

/// Perl's numeric coercion of a string: leading digits, or zero.
fn numeric(text: &str) -> usize {
    let digits: String = text.chars().take_while(|c| c.is_ascii_digit()).collect();
    digits.parse().unwrap_or(0)
}

/// Advance the read position past one MD segment that precedes an N.
///
/// A segment is either a run of digits (matching bases) or a mix of digits and mismatch
/// bases, possibly containing a `^` deletion.
fn advance(segment: &str) -> usize {
    if !segment.chars().any(|c| !c.is_ascii_digit()) {
        // A plain run of matches.
        return numeric(segment);
    }

    if segment.contains('^') {
        return advance_with_deletion(segment);
    }

    // Mismatch bases separated by match counts. Each mismatch base is one read base, hence
    // the `- 1` for the number of separators.
    let parts = perl_split(segment, |c| matches!(c, 'A' | 'T' | 'C' | 'G'));
    let sum: usize = parts.iter().map(|p| numeric(p)).sum();
    sum + parts.len().saturating_sub(1)
}

/// A segment containing a deletion: `^` followed by the deleted bases.
///
/// Deleted bases are present in the MD tag but absent from the read, so they do not advance
/// the read position, while mismatch bases do.
fn advance_with_deletion(segment: &str) -> usize {
    let mut pos = 0usize;
    let mut operation: Option<String> = None;

    for c in segment.chars() {
        let Some(current) = operation.as_mut() else {
            operation = Some(c.to_string());
            continue;
        };

        let current_is_number = current.chars().next().is_some_and(|f| f.is_ascii_digit());
        if current_is_number {
            if c.is_ascii_digit() {
                current.push(c);
            } else {
                pos += numeric(current);
                operation = Some(c.to_string());
            }
        } else if c.is_ascii_digit() {
            if !current.starts_with('^') {
                // A run of mismatch bases: each one is a read base.
                pos += current.chars().count();
            }
            operation = Some(c.to_string());
        } else {
            current.push(c);
        }
    }

    // The Perl requires the segment to end in a number and dies otherwise.
    pos + operation.as_deref().map(numeric).unwrap_or(0)
}

/// The first CIGAR operation this tool does not handle, if there is one.
///
/// Only M, I, D, S and N are accepted. `=` and `X` are a legal spelling of a match and a
/// mismatch, and aligners do emit them, but the Perl rejects them outright and a fixture pins
/// the message, so they are rejected here too rather than quietly supported.
pub fn unsupported_operation(cigar: &str) -> Option<char> {
    let ops = parse_cigar(cigar)?;
    ops.iter()
        .find(|op| !matches!(op.kind, b'M' | b'I' | b'D' | b'S' | b'N'))
        .map(|op| char::from(op.kind))
}

/// The read positions of every N in `md`, and their genomic positions.
///
/// `start` is the 1-based alignment position. Returned read positions index the read
/// sequence directly; returned genomic positions index the reference.
pub fn n_positions(md: &str, cigar: &str, start: usize) -> Option<Vec<(usize, usize)>> {
    let ops = parse_cigar(cigar)?;
    let segments: Vec<&str> = md.split('N').collect();
    if segments.len() < 2 {
        return Some(Vec::new());
    }

    let simple_match = ops.len() == 1 && ops[0].kind == b'M';
    let mut found = Vec::new();
    let mut pos = 0usize;
    let mut insertions_seen: Vec<usize> = Vec::new();
    let mut softclips_seen: Vec<usize> = Vec::new();

    // The last segment is whatever follows the final N and cannot itself precede one.
    for segment in &segments[..segments.len() - 1] {
        pos += advance(segment);

        let genomic = if simple_match {
            pos + start
        } else {
            pos = adjust_for_unrecorded_operations(
                pos,
                &ops,
                &mut insertions_seen,
                &mut softclips_seen,
            );
            genomic_position(pos, &ops, start)?
        };

        found.push((pos, genomic));
        // Step over the N itself.
        pos += 1;
    }

    Some(found)
}

/// Add the read bases that the MD tag never mentions.
///
/// Insertions and soft-clips consume read bases without appearing in the MD tag, so a
/// position derived from the MD tag alone is short by their length. Each one is counted at
/// most once, and only when it lies before the position being resolved.
fn adjust_for_unrecorded_operations(
    mut pos: usize,
    ops: &[Op],
    insertions_seen: &mut Vec<usize>,
    softclips_seen: &mut Vec<usize>,
) -> usize {
    let mut consumed = 0usize;

    for op in ops {
        match op.kind {
            b'M' => consumed += op.len,
            b'D' | b'N' => {}
            b'I' | b'S' => {
                if consumed >= pos {
                    // The operation is past the position being resolved.
                    break;
                }
                let seen = if op.kind == b'I' {
                    &mut *insertions_seen
                } else {
                    &mut *softclips_seen
                };
                if seen.contains(&consumed) {
                    continue;
                }
                seen.push(consumed);
                pos += op.len;
            }
            _ => {}
        }
    }

    pos
}

/// Walk the CIGAR to turn a read position into a genomic one.
fn genomic_position(pos: usize, ops: &[Op], start: usize) -> Option<usize> {
    let mut genomic = start;
    let mut remaining = pos as isize;

    for op in ops {
        match op.kind {
            b'M' => {
                if remaining <= op.len as isize {
                    return Some((genomic as isize + remaining) as usize);
                }
                genomic += op.len;
                remaining -= op.len as isize;
            }
            // Present in the read but not in the reference.
            b'I' | b'S' => remaining -= op.len as isize,
            // Present in the reference but not in the read.
            b'D' | b'N' => genomic += op.len,
            _ => return None,
        }
    }

    Some(genomic)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_cigar_splits_into_operations() {
        assert_eq!(
            parse_cigar("9M2I40M").unwrap(),
            vec![
                Op { kind: b'M', len: 9 },
                Op { kind: b'I', len: 2 },
                Op {
                    kind: b'M',
                    len: 40
                },
            ]
        );
        assert_eq!(
            parse_cigar("50M").unwrap(),
            vec![Op {
                kind: b'M',
                len: 50
            }]
        );
        assert!(parse_cigar("M50").is_none());
    }

    /// Perl drops trailing empty fields and the count is used in the arithmetic, so a Rust
    /// split would shift every position after a mismatch at the end of a segment.
    #[test]
    fn splitting_follows_perl_and_drops_trailing_empty_fields() {
        assert_eq!(perl_split("5A3", |c| c == 'A'), vec!["5", "3"]);
        assert_eq!(perl_split("5A", |c| c == 'A'), vec!["5"]);
        assert_eq!(perl_split("A5", |c| c == 'A'), vec!["", "5"]);
    }

    #[test]
    fn a_plain_match_run_advances_by_its_length() {
        assert_eq!(advance("24"), 24);
        assert_eq!(advance("0"), 0);
    }

    /// `10A5` is ten matches, one mismatch base, five matches: sixteen read bases.
    #[test]
    fn a_mismatch_base_counts_as_one_read_base() {
        assert_eq!(advance("10A5"), 16);
        assert_eq!(advance("10A5T2"), 19);
    }

    /// Deleted bases are in the MD tag but not in the read.
    #[test]
    fn deleted_bases_do_not_advance_the_read_position() {
        assert_eq!(advance("33^ACTCGA11"), 44);
    }

    /// The simple single-match path: MD:Z:24N25 on a 50M read starting at 6.
    #[test]
    fn a_single_match_read_resolves_its_n_directly() {
        let found = n_positions("24N25", "50M", 6).unwrap();
        assert_eq!(found, vec![(24, 30)]);
    }

    #[test]
    fn two_ns_in_one_read_are_both_found() {
        let found = n_positions("10N29N9", "50M", 20).unwrap();
        assert_eq!(found.len(), 2);
        assert_eq!(found[0], (10, 30));
        // The second N is 29 further on, plus one for the first N itself.
        assert_eq!(found[1], (40, 60));
    }

    /// An insertion consumes read bases the MD tag never mentions, so the read position has
    /// to be pushed past it while the genomic position is not.
    #[test]
    fn an_insertion_shifts_the_read_position_but_not_the_genomic_one() {
        // 9M2I40M, MD:Z:44N4: the N is at read position 46 and genomic 44 past the start.
        let found = n_positions("44N4", "9M2I40M", 31124195).unwrap();
        let (read_pos, genomic) = found[0];
        assert_eq!(read_pos, 46);
        assert_eq!(genomic, 31124195 + 44);
    }

    /// A deletion is the mirror image: the genomic position advances, the read position does
    /// not.
    #[test]
    fn a_deletion_shifts_the_genomic_position_but_not_the_read_one() {
        // 33M6D18M, MD:Z:33^ACTCGA11N6.
        let found = n_positions("33^ACTCGA11N6", "33M6D18M", 87082089).unwrap();
        let (read_pos, genomic) = found[0];
        assert_eq!(read_pos, 44);
        assert_eq!(genomic, 87082089 + 44 + 6);
    }

    #[test]
    fn an_md_tag_with_no_n_yields_nothing() {
        assert_eq!(n_positions("50", "50M", 1).unwrap(), vec![]);
    }
}
