# Distribution and documentation changes

Scientific results: corrected Yanayacu reviewer edition, 8 October 2026. Computational distribution: 26 November 2026.

## Repository documentation revision

This revision shortens six documentation files and corrects the usage example in the `figure_code/build_alpha.py` module docstring to use the public figure-rendering entry point. Executable statements, source data, numerical results and figures are unchanged. No scientific analyses were rerun for this editorial revision.

The ZIP retains its original documentation and script bytes; the GitHub branch contains the later editorial revision. The archive checksum remains in `DOWNLOAD_MANIFEST.json`, and `SHA256_MANIFEST.tsv` describes that archived distribution.

## Computational distribution

- Includes recorded and prepared metadata, retained ASV and genus-lineage matrices, statistical tables, publication figures and native execution records. FASTQ reads are linked through the deposited sequence record and sample manifest.
- Reproduces ten association tables byte for byte. Figure-code portability changes preserve source values, and the distributed publication figures retain their original bytes.
- Relocates table-index links to the standalone package without changing workbook values.
- Retains script identifiers, execution dates, input preparation records and source hashes. Earlier complete deliveries, backup ZIPs and confidential reviewer material are excluded.

A complete rerun from FASTQ was not performed. The reproduction guides distinguish the executed numerical checks from retained native objects and stages requiring the original environments and reads.
