# Build the suite, then ship one static-ish binary and its three classic names.
#
# The image carries no Perl and no samtools: that is the point of the rewrite, and an image
# that quietly installed them would hide whether the claim is true.

FROM rust:1.89-slim AS build
WORKDIR /src
COPY . /src
RUN cargo build --release --locked

FROM debian:bookworm-slim
LABEL org.opencontainers.image.title="SNPsplit"
LABEL org.opencontainers.image.description="Allele-specific alignment sorting"
LABEL org.opencontainers.image.source="https://github.com/FelixKrueger/SNPsplit"
LABEL org.opencontainers.image.licenses="GPL-3.0-only"

# ca-certificates is needed by --download and by nothing else. Without it the tools still
# work; only fetching references would fail, and it would fail with a TLS error rather than
# something useful.
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates \
 && rm -rf /var/lib/apt/lists/*

COPY --from=build /src/target/release/snpsplit /usr/local/bin/snpsplit
# The classic names, so existing pipelines are drop-in. Each is a link to the one binary,
# which dispatches on argv[0].
RUN ln -s /usr/local/bin/snpsplit /usr/local/bin/SNPsplit \
 && ln -s /usr/local/bin/snpsplit /usr/local/bin/tag2sort \
 && ln -s /usr/local/bin/snpsplit /usr/local/bin/SNPsplit_genome_preparation

ENTRYPOINT ["/usr/local/bin/snpsplit"]
