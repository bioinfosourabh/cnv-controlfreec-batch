# Control-FREEC 11.6b with bedtools, needed for BAF on targeted (exome/panel)
# data: FREEC runs `bedtools intersect` to restrict SNP positions to the
# capture regions before building pileups.
#
#   docker build -t freec-bedtools:11.6b containers/freec-bedtools
#   nextflow run . ... --freec_container freec-bedtools:11.6b
FROM quay.io/biocontainers/control-freec:11.6b--hde5307d_3

ADD --chmod=755 https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static /usr/local/bin/bedtools

RUN bedtools --version
