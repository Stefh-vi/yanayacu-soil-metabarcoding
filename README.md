# Yanayacu soil metabarcoding: current analysis and source data

This package contains the current 16S and ITS soil-metabarcoding analysis for nine soil units from Jaén and San Ignacio, Peru. The scientific snapshot corresponds to the corrected reviewer edition of 8 October 2026. This standalone distribution was assembled on 26 November 2026 without copying earlier delivery folders or backup archives.

For the complete source data, native statistical objects and publication figures, download and extract `YANAYACU_CURRENT_CODE_DATA_20261126.zip` from the [current complete distribution](https://drive.google.com/drive/folders/1mZ3jRF7mkHfaLVXjxxBNuF8gpr4TbRUf). The folder contains the ZIP, its SHA-256 checksum and verification receipts. GitHub exposes the executable source and documentation for inspection. Run the commands below from the fully extracted distribution, not from a source-only checkout. The archive's checksum and file manifest identify the exact delivered contents. A GitHub release asset and a Zenodo DOI have not been published.

Both markers describe the same nine analytical units: B1–B3 (low degradation), M1–M3 (moderate degradation), and A1–A3 (high degradation). J and SI retain their original locality meanings. The original fertile-soil unit P1 is outside the comparative cohort.

## Start here

1. Read [the analysis reproduction guide](docs/ANALYSIS_REPRODUCTION_EN.md).
2. Read [the raw data and sequence workflow guide](docs/RAW_DATA_AND_SEQUENCE_WORKFLOW_EN.md) for original read access, sample mapping and native workflow records.
3. Read [the figure reproduction guide](docs/FIGURE_REPRODUCTION_EN.md) for the current figure builders and preserved native figures.
4. Consult [provenance](PROVENANCE.md), the source manifests and validation reports for the scope of each check.

## Contents

| Directory | Contents |
|---|---|
| `data/metadata/` | Original recorded and analysis-prepared sample metadata and identity mapping |
| `data/sequence_inputs/` | Retained complete ASV matrices and associated sequence-analysis source data |
| `data/raw_manifest/` | Original FASTQ identity and public sequence-accession/download mapping |
| `data/analysis_inputs/` | Complete current statistical input matrices and frozen analysis specification |
| `data/figure_inputs/` | Numerical and native source inputs required by the figure builders |
| `analysis/` | Current statistical code, independent verification and reference output hashes |
| `sequence_workflow/` | Source scripts, commands, parameters and environments for the retained native workflow |
| `figure_code/` | Portable current figure rendering entry point |
| `tables/` | Current publication workbook and machine-readable result tables |
| `figures/` | Current publication figures and their captions |
| `validation/` | Execution-specific numerical, graphical, raw-data and file-identity receipts |
| `docs/` | Reproduction instructions and current figure/table legends |

## Verify the downloaded files

From the package root, with Python 3.10 or later:

```sh
python scripts/verify_release.py --report ../YANAYACU_FILE_VERIFICATION.json
```

The report must be outside the released package. Work in a copy and keep the original download immutable.

## Reproduce the current numerical analyses

Use a separate environment. The statistical replay uses Python 3.12.14 and NumPy 2.3.5:

```sh
python -m pip install -r analysis/requirements.txt
python analysis/replay.py --output ../YANAYACU_ANALYSIS_REPLAY
python analysis/verify_metadata.py --report ../YANAYACU_ANALYSIS_REPLAY/validation/METADATA_PREPARATION_QA.json
```

Choose a new output directory. The verified public replay reproduces ten current association tables byte for byte. It independently checks four alpha-category models and 28 environmental category comparisons. Retained ANCOM-BC2 and dispersion arithmetic checks are not native model refits.

## Reproduce the current figures

Follow `docs/FIGURE_REPRODUCTION_EN.md` for the rendering environment, then run:

```sh
python figure_code/reproduce_figures.py --output ../YANAYACU_FIGURE_REPLAY
```

Current publication files are supplied in `figures/`. The replay distinguishes figures regenerated from numerical tables, native figures with revised annotations, and current native figures retained from their source. It does not refit scientific models. Rendering-library versions may change font metrics or PDF bytes; the validation records separately compare numerical labels, page geometry and rendered appearance.

## Original reads and the master data

The deposited sequence records are associated with NCBI BioProject [PRJNA1455726](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1455726). Raw FASTQ files remain in the sequence archive and are referenced by the supplied accession and checksum manifest. They are not embedded in this code package. The original deposition contains P1 in addition to the nine analytical units; use the mapping and inclusion fields rather than treating every project run as an analysis sample.

The original recorded metadata and current prepared metadata are preserved separately. Complete retained ASV and genus-lineage matrices are processed source data, not raw reads. All input data required by the standalone statistical replay are included; it does not require a workstation account or a previous delivery directory.

Read the raw-data guide before downloading or remapping reads. Publicly recompressed FASTQ files can differ in compressed-byte checksum from vendor files. Discrepancies in public experiment titles are documented and must not be silently interpreted as corrected repository metadata.

## Scope and interpretation

The November packaging does not constitute a new complete analysis from FASTQ. Native bioinformatic stages require the actual matching raw reads, reference resources, recorded environments and path configuration. Their executed source records are retained. The numerical replay, figure rendering and raw-data lookup are distinct operations with separate validation reports.

The study has nine shared soil units, not eighteen independent marker samples or one hundred biological ITS replicates. Low degradation is absent from San Ignacio. Source files do not establish sampling independence, complete chemical assay documentation or causal relationships. Original and prepared particle-size values retain their documented distinction.

No earlier complete deliverables, confidential reviewer correspondence or backup ZIPs are part of this distribution. Required source inputs and scripts that predate the scientific snapshot remain traceable through their source manifests.

