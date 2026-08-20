//! Scoring a masked position in bisulfite data.
//!
//! Bisulfite converts unmethylated C to T on the strand that was read, so a C in the genome
//! may legitimately appear as either C or T. That makes some SNPs unusable on some strands:
//! a C>T SNP is indistinguishable from an unconverted versus converted C on the top strand,
//! and a G>A SNP is the same problem on the bottom strand.

use super::snps::Snp;

/// Which of the four bisulfite strands a read came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// Original top.
    Ot,
    /// Complementary to original top.
    CtOt,
    /// Original bottom.
    Ob,
    /// Complementary to original bottom.
    CtOb,
}

impl Strand {
    /// Derive the strand from Bismark's `XR:Z:` and `XG:Z:` tags.
    pub fn from_conversions(read: &str, genome: &str) -> Option<Self> {
        match (read, genome) {
            ("CT", "CT") => Some(Strand::Ot),
            ("GA", "CT") => Some(Strand::CtOt),
            ("GA", "GA") => Some(Strand::CtOb),
            ("CT", "GA") => Some(Strand::Ob),
            _ => None,
        }
    }

    fn is_top(self) -> bool {
        matches!(self, Strand::Ot | Strand::CtOt)
    }
}

/// What scoring one position concluded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Call {
    Genome1,
    Genome2,
    /// The base matched neither allele.
    Neither,
    /// The position could not be used because bisulfite conversion makes the two alleles
    /// indistinguishable on this strand.
    UnusableCtSnp,
}

/// Score one masked position for a bisulfite read.
///
/// The four branches are the four ways a SNP can involve a cytosine: reference C,
/// alternative C, reference G, alternative G. A SNP touching neither is scored as it would be
/// in a normal library.
pub fn score(snp: Snp, read_base: u8, strand: Strand) -> Call {
    let (reference, alternative) = (snp.reference, snp.alternative);

    // Reference is C: on the top strand a C in the read may have been a converted T.
    if reference == b'C' {
        if alternative == b'T' && strand.is_top() {
            return Call::UnusableCtSnp;
        }
        return if strand.is_top() {
            match read_base {
                b'C' | b'T' => Call::Genome1,
                b if b == alternative => Call::Genome2,
                _ => Call::Neither,
            }
        } else if read_base == reference {
            Call::Genome1
        } else if alternative == b'G' {
            // The bottom strand's own conversion: G may appear as A.
            match read_base {
                b'G' | b'A' => Call::Genome2,
                _ => Call::Neither,
            }
        } else if read_base == alternative {
            Call::Genome2
        } else {
            Call::Neither
        };
    }

    // Alternative is C: the mirror image, with the genomes swapped.
    if alternative == b'C' {
        if reference == b'T' && strand.is_top() {
            return Call::UnusableCtSnp;
        }
        return if strand.is_top() {
            match read_base {
                b'C' | b'T' => Call::Genome2,
                b if b == reference => Call::Genome1,
                _ => Call::Neither,
            }
        } else if read_base == alternative {
            Call::Genome2
        } else if reference == b'G' {
            match read_base {
                b'G' | b'A' => Call::Genome1,
                _ => Call::Neither,
            }
        } else if read_base == reference {
            Call::Genome1
        } else {
            Call::Neither
        };
    }

    // Reference is G: the cytosine is on the opposite strand, so the bottom strand converts.
    if reference == b'G' {
        if alternative == b'A' && !strand.is_top() {
            return Call::UnusableCtSnp;
        }
        return if strand.is_top() {
            if read_base == reference {
                Call::Genome1
            } else if read_base == alternative {
                Call::Genome2
            } else {
                Call::Neither
            }
        } else {
            match read_base {
                b'G' | b'A' => Call::Genome1,
                b if b == alternative => Call::Genome2,
                _ => Call::Neither,
            }
        };
    }

    // Alternative is G.
    if alternative == b'G' {
        if reference == b'A' && !strand.is_top() {
            return Call::UnusableCtSnp;
        }
        return if strand.is_top() {
            if read_base == alternative {
                Call::Genome2
            } else if read_base == reference {
                Call::Genome1
            } else {
                Call::Neither
            }
        } else {
            match read_base {
                b'G' | b'A' => Call::Genome2,
                b if b == reference => Call::Genome1,
                _ => Call::Neither,
            }
        };
    }

    // The SNP touches no cytosine on either strand.
    if read_base == reference {
        Call::Genome1
    } else if read_base == alternative {
        Call::Genome2
    } else {
        Call::Neither
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn snp(reference: u8, alternative: u8) -> Snp {
        Snp {
            reference,
            alternative,
        }
    }

    #[test]
    fn the_strand_comes_from_the_two_bismark_tags() {
        assert_eq!(Strand::from_conversions("CT", "CT"), Some(Strand::Ot));
        assert_eq!(Strand::from_conversions("GA", "CT"), Some(Strand::CtOt));
        assert_eq!(Strand::from_conversions("GA", "GA"), Some(Strand::CtOb));
        assert_eq!(Strand::from_conversions("CT", "GA"), Some(Strand::Ob));
        assert_eq!(Strand::from_conversions("CT", "XX"), None);
    }

    /// The central case: on the top strand a C>T SNP is indistinguishable from an
    /// unconverted versus converted cytosine, so the position tells you nothing.
    #[test]
    fn a_c_to_t_snp_is_unusable_on_the_top_strand() {
        assert_eq!(
            score(snp(b'C', b'T'), b'T', Strand::Ot),
            Call::UnusableCtSnp
        );
        assert_eq!(
            score(snp(b'C', b'T'), b'T', Strand::CtOt),
            Call::UnusableCtSnp
        );
    }

    /// The same SNP is perfectly usable on the bottom strand, where no conversion applies.
    #[test]
    fn a_c_to_t_snp_is_usable_on_the_bottom_strand() {
        assert_eq!(score(snp(b'C', b'T'), b'C', Strand::Ob), Call::Genome1);
        assert_eq!(score(snp(b'C', b'T'), b'T', Strand::Ob), Call::Genome2);
    }

    #[test]
    fn a_g_to_a_snp_is_the_mirror_image() {
        assert_eq!(
            score(snp(b'G', b'A'), b'A', Strand::Ob),
            Call::UnusableCtSnp
        );
        assert_eq!(score(snp(b'G', b'A'), b'G', Strand::Ot), Call::Genome1);
        assert_eq!(score(snp(b'G', b'A'), b'A', Strand::Ot), Call::Genome2);
    }

    /// A C that is not part of a C>T SNP still converts, so both C and T mean genome 1.
    #[test]
    fn a_reference_c_accepts_its_converted_form_on_the_top_strand() {
        assert_eq!(score(snp(b'C', b'A'), b'C', Strand::Ot), Call::Genome1);
        assert_eq!(score(snp(b'C', b'A'), b'T', Strand::Ot), Call::Genome1);
        assert_eq!(score(snp(b'C', b'A'), b'A', Strand::Ot), Call::Genome2);
        assert_eq!(score(snp(b'C', b'A'), b'G', Strand::Ot), Call::Neither);
    }

    #[test]
    fn a_snp_touching_no_cytosine_is_scored_normally() {
        assert_eq!(score(snp(b'A', b'T'), b'A', Strand::Ot), Call::Genome1);
        assert_eq!(score(snp(b'A', b'T'), b'T', Strand::Ob), Call::Genome2);
        assert_eq!(score(snp(b'A', b'T'), b'G', Strand::Ot), Call::Neither);
    }
}
