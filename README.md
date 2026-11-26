# Yanayacu soil metabarcoding

16S and ITS analysis of nine soil units from Jaén and San Ignacio, Peru. The results correspond to the corrected reviewer edition of 8 October 2026; the computational distribution was assembled on 26 November 2026.

## Download the data and code

Download and extract `YANAYACU_CURRENT_CODE_DATA_20261126.zip` from the [complete distribution](https://drive.google.com/drive/folders/1mZ3jRF7mkHfaLVXjxxBNuF8gpr4TbRUf). It contains the analysis inputs, native statistical objects, tables, figures and validation reports. The archive size and checksum are recorded in [DOWNLOAD_MANIFEST.json](DOWNLOAD_MANIFEST.json).

Run the commands below from the extracted archive. The GitHub checkout provides source code and documentation but omits some inputs required for execution. No GitHub release asset or Zenodo DOI has been published.

GitHub includes a later editorial revision to the guides and an alpha-renderer usage example. The ZIP retains its original files; `SHA256_MANIFEST.tsv` applies to that archive. See [RELEASE_NOTES.md](RELEASE_NOTES.md) for the revision scope.

Both markers characterize the same units: B1-B3 (low degradation), M1-M3 (moderate) and A1-A3 (high). J denotes Jaén and SI denotes San Ignacio. P1, the fertile-soil reference, is excluded from the comparative cohort.

## Package contents

| Directory | Contents |
|---|---|
| `data/metadata/` | Recorded measurements, prepared analysis values and sample identifiers |
| `data/sequence_inputs/` | Retained ASV counts and taxonomy |
| `data/raw_manifest/` | FASTQ filenames, accessions and checksums |
| `data/analysis_inputs/` | Statistical input matrices and model specification |
| `data/figure_inputs/` | Plotting data and native figure inputs |
| `analysis/` | Statistical code, independent checks and reference hashes |
| `sequence_workflow/` | Executed source records, parameters and environments |
| `figure_code/` | Figure rendering scripts |
| `tables/` | Publication workbook and machine-readable tables |
| `figures/` | Publication figures and captions |
| `validation/` | Numerical, figure and file-identity checks |
| `docs/` | Reproduction guides and figure/table legends |

## Verify the download

From the package root, with Python 3.10 or later:

```sh
python scripts/verify_release.py --report ../YANAYACU_FILE_VERIFICATION.json
```

Keep the verification report and subsequent outputs outside the extracted distribution.

## Reproduce the statistical analyses

Use a separate environment. The verified statistical replay used Python 3.12.14 and NumPy 2.3.5. Choose a new output directory:

```sh
python -m pip install -r analysis/requirements.txt
python analysis/replay.py --output ../YANAYACU_ANALYSIS_REPLAY
python analysis/verify_metadata.py --report ../YANAYACU_ANALYSIS_REPLAY/validation/METADATA_PREPARATION_QA.json
```

This replay regenerated ten association tables byte for byte and independently checked four alpha-category models and 28 environmental category comparisons. The ANCOM-BC2 and dispersion checks verify retained table arithmetic; they do not refit the native models. See the [statistical reproduction guide](docs/ANALYSIS_REPRODUCTION_EN.md) for model specifications, correction families and validation results.

## Reproduce the figures

Install the environment specified in the [figure reproduction guide](docs/FIGURE_REPRODUCTION_EN.md), then run:

```sh
python figure_code/reproduce_figures.py --output ../YANAYACU_FIGURE_REPLAY
```

The script renders plots from numerical tables, updates two native figure annotations and copies the remaining publication figures. It checks numerical labels, page geometry and appearance. Font and library differences can change rebuilt PDF bytes. Statistical models are not refitted by this command.

## Raw reads and metadata

Raw FASTQ files are available through [BioProject PRJNA1455726](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1455726) and the supplied accession manifest; they are external to the ZIP. The deposition includes P1, so use the inclusion fields to select the nine-unit cohort. Public and vendor FASTQ representations have separate checksums. Discrepancies between public experiment titles and library names remain unresolved in the deposited metadata; the [raw-data guide](docs/RAW_DATA_AND_SEQUENCE_WORKFLOW_EN.md) explains the sample mapping and identity checks.

The archive includes all inputs for the standalone statistical replay. The ASV and genus-lineage matrices contain processed counts. Recorded measurements and prepared metadata are preserved separately, including the A1 particle-size normalization from 102% to 100%; that transformation is not a laboratory correction.

## Scope and provenance

The November checks did not rerun the complete FASTQ-to-figure workflow. Repeating the recorded sequence stages requires the matching reads, reference resources, original environments and path configuration. Their source records are retained in `sequence_workflow/`.

The study comprises nine shared soil units. The two markers and the ITS technical subsamples do not increase the biological sample size. Low degradation occurs only in Jaén; field independence and complete chemical assay documentation remain unverified. These observational data do not establish causal effects.

[PROVENANCE.md](PROVENANCE.md) records input preparation, source history and computational assistance. The distribution retains the source scripts and hashes required to trace the results, while excluding earlier complete deliveries, backup archives and confidential reviewer correspondence.
