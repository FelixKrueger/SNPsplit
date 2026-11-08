# Build the suite, then ship one binary and its three classic names.
#
# The runtime image is distroless: no shell, no package manager, and no Perl. That last one
# is not decoration. Every Debian base image includes perl-base as an essential package, so
# an image built on debian:*-slim ships a Perl interpreter no matter what this file says, and
# an image making a "no Perl" claim on that base would be making it falsely.

# Built on bookworm rather than the default slim tag so the binary links against the same
# glibc the runtime image has. Building on a newer base produced a binary that failed at
# startup with `GLIBC_2.39 not found`, which is the kind of thing only running the image
# catches.
FROM rust:1.89-bookworm AS build
WORKDIR /src
COPY . /src
RUN cargo build --release --locked

# The classic names are made here because the runtime image has no shell to make them in.
RUN mkdir -p /out/bin \
 && cp target/release/snpsplit /out/bin/snpsplit \
 && ln -s snpsplit /out/bin/SNPsplit \
 && ln -s snpsplit /out/bin/tag2sort \
 && ln -s snpsplit /out/bin/SNPsplit_genome_preparation

# cc rather than static: the binary links glibc. This base carries ca-certificates, which
# --download needs and nothing else does.
FROM gcr.io/distroless/cc-debian12
LABEL org.opencontainers.image.title="SNPsplit"
LABEL org.opencontainers.image.description="Allele-specific alignment sorting"
LABEL org.opencontainers.image.source="https://github.com/FelixKrueger/SNPsplit"
LABEL org.opencontainers.image.licenses="GPL-3.0-only"

COPY --from=build /out/bin /usr/local/bin
ENV PATH="/usr/local/bin:${PATH}"

ENTRYPOINT ["/usr/local/bin/snpsplit"]
