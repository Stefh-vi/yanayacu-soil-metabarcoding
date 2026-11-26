# Raw data, master metadata and sequence provenance

This release packages the inputs and active processing records underlying the current nine-sample analysis. It does not contain earlier delivery editions or their result archives. Earlier execution dates in source records identify how the retained results were produced; they are not alternative current analyses.

## What is included

| Layer | Location | Scope |
|---|---|---|
| Current metadata | `data/metadata/soil_metadata_analysis.tsv` | Exact current S16 table: nine sampling units |
| Metadata preparation record | `data/metadata/soil_metadata_source_with_preparation.tsv` | Source values, prepared texture values, original texture values and preparation label |
| Identity crosswalk | `data/metadata/sample_identity.tsv` | B1–B3, M1–M3, A1–A3; vendor identifiers; J = Jaen and SI = San Ignacio |
| Retained ASV counts and taxonomy | `data/sequence_inputs/` | Exact upstream matrices used for the current analysis, with taxonomy and 16S filtering decisions |
| Public raw-data crosswalk | `data/raw_manifest/` | Runs, original vendor filenames/hashes, ENA filenames/MD5/size and canonical Drive links |
| Sequence processing evidence | `sequence_workflow/recorded_sources/`, `parameters/`, `environments/` | Active source code, recorded arguments and installed environments |
| Native statistical analysis | `sequence_workflow/native_analysis/` | Final ANCOM-BC2 R objects, frozen input counts, R code, environment locks, beta inputs and native NMDS source/objects |
| Quality reports | `sequence_workflow/qc_reports/` | MultiQC after-trimming reports for both markers; download and transformation receipt under `validation/raw/` |

The source metadata file is a preparation record, not an original laboratory certificate. For A1, the original sand/silt/clay values were 68.82/17.18/16.00, summing to 102. The prepared values close the composition to 100; the original values and preparation record remain available. This is not evidence of a laboratory correction. Total nitrogen, organic matter and total carbon remain in the metadata. Laboratory units, analytical methods and the relationship between reported carbon and organic matter must not be invented.

## Public raw data and identity checks

The official [ENA file-report service](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1455726&result=read_run&fields=run_accession,library_name,experiment_title,read_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv) returned 20 runs for [BioProject PRJNA1455726](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1455726) on 26 November 2026. Eighteen runs represent nine retained sampling units with two markers. The two P1660/P1 runs are recorded as excluded reference material and do not enter the current nine-sample analyses.

The 36 ENA FASTQ files for the analysis have a total declared compressed size of **1,686,854,967 bytes**. The canonical Drive folder contains 36 vendor-named FASTQ files with a total declared size of **1,601,546,685 bytes**. Deposited read counts agree with all 18 current raw read-pair counts: **4,722,910 pairs for 16S** and **6,529,501 pairs for ITS**. These agreements verify the recorded library/count mapping, not nucleotide-level or compressed-byte identity. No FASTQ files were downloaded or rehashed during this packaging operation. Canonical vendor SHA-256 values come from the archived source manifest; ENA MD5 values describe the public archive files. They are different checksum domains.

Eighteen of the 20 public experiment titles have a sample prefix inconsistent with their library name. The library names and canonical vendor/sample crosswalk provide the mapping retained here; the title discrepancy is explicitly recorded in `public_accession_map.tsv`. It requires depositor reconciliation and must not be silently treated as resolved. Exact public metadata and the HTTP retrieval receipt are retained in `validation/raw/`.

The author's canonical source folder is [RAW data on Drive](https://drive.google.com/drive/folders/10mN-aWoFK-vfuwB1vJRj7dN0vxjhjoMV). Drive access conditions may differ from the public ENA record. The raw bytes are not duplicated inside this code release.

To inspect the included download manifest without network transfer:

```sh
python sequence_workflow/download_public_fastq.py
```

To download the mapped public FASTQ files into a **new** directory and verify each public file's size and MD5:

```sh
python sequence_workflow/download_public_fastq.py --download --output /path/to/new/yanayacu_fastq
```

`--marker 16S` or `--marker ITS` restricts the transfer. `--sample B1` can be repeated for a subset. This command has passed manifest-only preflight; a download was not run for this release. Its download receipt will document new SHA-256 values separately from vendor checksums.

## Sequence processing and reproducibility boundary

The active sequence algorithms are included as `.py.txt` records. Host-specific project and environment prefixes were replaced by named placeholders; original and packaged SHA-256 values are listed in `validation/raw/SOURCE_FILE_INVENTORY.json`. These files are reference source records, not a validated turnkey launcher. They include the real processing decisions and command arguments rather than a newly invented workflow.

The retained 16S settings include Cutadapt 5.1, operational primer candidates `ACTCCTACGGGAGGCAGCAG` and `GGACTACHVGGGTWTCTAAT`, 270/210 forward/reverse truncation, expected-error limits 2/2, 20-base minimum overlap, zero merge mismatches and consensus chimera removal. The ITS branch uses operational ITS3/ITS4-compatible cores `GCATCGATGAAGAACGCAGC` and `TCCTCCGCTTATTGATATGC`, preserves variable length with fixed truncation disabled, and records its own filtering, overlap and chimera settings. These operational oligonucleotides are not claimed to be vendor-certified primer sequences.

The final ecology datasets contain **27,834 retained 16S ASVs / 1,752,112 sequences** and **8,927 retained ITS ASVs / 5,388,961 sequences**, each across nine samples. Final alpha/rarefaction depths are **174,661 for 16S** and **450,343 for ITS**. The 16S alpha estimator is one rarefaction realization; the ITS estimator is the mean of 100 technical subsamples. Technical subsampling intervals are not biological confidence intervals. The parameter record also contains pre-taxonomy technical-stage settings; these are identified by scope and must not replace final ecological settings.

A complete raw-to-figure rerun was **not** performed for this packaging. Restoring the original execution environments, reference classifier/database assets, path mappings and raw inputs is still required to replay the recorded sequence stages. The exact original 16S rarefied ASV realization was not found in the selected local source archive; its retained distance matrix, alpha outcomes and rarefaction summary are included. A new draw from the unrarefied counts must not be presented as the original draw. No new Nextflow execution or end-to-end validation is claimed.

## Native ANCOM-BC2 replay

`sequence_workflow/native_analysis/ancombc2/` contains the final extended native fit only: 356 eligible 16S full-lineage features and 179 eligible ITS full-lineage features, formula `locality + zone`, marker-specific global Holm families, and native pseudocount sensitivity. The base frozen configuration and numerical-extension addendum jointly define the effective configuration: EM iteration ceiling **10,000**, tolerance **1e-5**. The addendum records why the iteration budget changed; initial unconverged result files are not included.

The exact R execution source, frozen input hashes, native RDS objects, coefficient/global result tables, convergence diagnostics and R session information are included. Native global W/raw-P/Holm-P values match all **535 current S17 rows exactly**. The package preflight verifies **14 source/input hashes and 228 pinned Linux packages**. These checks validate the packaged lineage and byte consistency; they are not a fresh fit.

```sh
python sequence_workflow/native_analysis/ancombc2/replay_from_package.py --check-only
```

On a compatible Linux host, restore the recorded Conda environment using the explicit lock, then fit into a new output directory:

```sh
conda create --prefix /path/to/new/ancombc2_env --file sequence_workflow/native_analysis/ancombc2/executed_run/CONDA_EXPLICIT_WITH_MD5.txt
python sequence_workflow/native_analysis/ancombc2/replay_from_package.py --rscript /path/to/new/ancombc2_env/bin/Rscript --output /path/to/new/ancombc2_results
```

The launcher requires **R 4.5.3 / ANCOMBC 2.14.0** and checks hashes before fitting. This native refit was not executed in the publication packaging environment. Native source status fields may have broader all-coefficient validity flags; current S17/S18 explicitly distinguish the primary degradation endpoints from auxiliary coefficient validity. The current reporting tables remain the interpretation reference.

## Native beta, dispersion and NMDS

`native_analysis/ecology/inputs/` includes the exact retained 16S Bray/Lingoes and ITS Bray-proportion/Hellinger distance matrices, 4,319 nonidentity within-locality permutations, current category/dispersion reference tables, and native rarefaction summaries. The new `replay_beta.R` entry point refits three conditional beta models and three native PERMDISP diagnostics from these retained matrices. It checks native raw results against current tables and joins current alpha P values to enforce the final **16S Holm family of three endpoints and ITS Holm family of four**. The four alpha tests are independently refitted by the main numerical replay; this beta wrapper imports those current alpha reference values rather than refitting alpha.

```sh
Rscript --vanilla sequence_workflow/native_analysis/ecology/replay_beta.R /path/to/release /path/to/new/beta_replay
Rscript --vanilla sequence_workflow/native_analysis/ecology/run_nmds.R sequence_workflow/native_analysis/ecology/inputs /path/to/new/nmds_replay
```

Both require vegan 2.6.8 and its recorded dependencies. NMDS uses the exact retained source algorithm, fixed seeds and native source objects. The R entry points are supplied for replay but were **not executed or R-parser validated** during this packaging because the required R environment was unavailable locally. Archived beta source files remain `.R.txt` records: the original ITS source used a two-distance multiplicity family, which is superseded by the current four-endpoint family and must not be reused as current inference. Within-locality permutations do not make PERMDISP a locality-adjusted dispersion model.

## Quality report integrity

The included MultiQC HTML reports were obtained from the authorized Drive source. Each connector-delivered response was three bytes longer than Drive's declared file size. Their received SHA-256 values and byte counts are recorded; remote byte identity is therefore **not asserted**. Only the workstation project prefix in each exported report was replaced by `${WORKFLOW_ROOT}` for portability. Embedded QC values, figures and software notices were retained. The HTML document structure was checked; these copies were not rendered in a browser in this packaging pass. Original native reports and larger QC/QIIME archives remain linked in the canonical raw source folder.

## What the current checks establish

The raw library mapping, current read counts, metadata input hashes, retained ASV matrix totals, current ANCOM numerical outputs and packaged source locks are verified. They do not establish missing laboratory units/methods, provider primer certification, nucleotide identity of two archive representations, or causal effects of degradation. There are nine sampling units, low degradation is absent from SI, and locality/vegetation structure limits inference. These are properties of the evidence, not problems that can be removed by packaging.

See `validation/raw/` for machine-readable checks, `ANALYSIS_REPRODUCTION_EN.md` for numerical refits, `FIGURE_REPRODUCTION_EN.md` for current figure reproduction, and `sequence_workflow/THIRD_PARTY_NOTICES.md` for native software attribution.
