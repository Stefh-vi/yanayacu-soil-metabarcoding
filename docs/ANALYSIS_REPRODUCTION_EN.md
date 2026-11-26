# Reproduce the current statistical analyses

This package contains the numerical specification, current source data, analysis code and current reference outputs associated with the corrected Yanayacu reviewer edition. It does not require any previous delivery folder. The November packaging changes paths and removes comparison-only material; it does not change the released statistical estimates or conclusions.

## Quick start

Python 3.12.14 and NumPy 2.3.5 were used for the verified replay. From the extracted package directory:

```sh
python -m pip install -r analysis/requirements.txt
python analysis/replay.py --output ../YANAYACU_REPLAY_NEW
python analysis/verify_metadata.py --report ../YANAYACU_REPLAY_NEW/validation/METADATA_PREPARATION_QA.json
```

Choose an output directory that does not yet exist. The replay refuses to overwrite earlier outputs and writes outside the released package. It requires no workstation account, absolute source path, network service or raw FASTQ download. Installation is only needed when NumPy is not already available at the recorded version.

The replay executes `analysis/run_associations.py`, `analysis/verify_independently.py` and `analysis/verify_retained_statistics.py`. Current output tables are written under the new directory's `tables/`; numerical reports and category refits are written under `validation/`.

## Source data included

`data/analysis_inputs/` contains:

- `S16_Study_metadata.tsv`: all nine analytical units and the 28 recorded numerical soil variables, locality and vegetation. Original units that could not be established remain explicitly unconfirmed.
- `metadata_preparation_record.tsv`: the corresponding preparation record, including long/short sample ID mapping and original versus analytical particle-size values. This is provenance for the current input, not an obsolete analysis version.
- The complete retained genus-rank lineage count matrices for both markers: 1,143 prokaryotic lineages and 1,185 fungal lineages. These are processed counts; they are not raw sequencing reads and are not the full ASV tables.
- Current genus selection and per-sample data; 18 marker/sample alpha estimates; seven current alpha/beta category records; the 28-variable environmental comparison table; retained ANCOM-BC2 and dispersion tables.
- `DESIGN.json`: the current model and correction families, input hashes and sample order. The preserved October specification timestamp is not a claim that the study was preregistered. The November packaging timestamp and source specification hash are recorded separately.

All copied numerical input files and the ten current reference output tables retain the exact source bytes. `analysis/SOURCE_MANIFEST.json` records their source-edition members and SHA-256 hashes without personal filesystem paths. A single source version is used throughout.

## Statistical scope

| Component | Reproduction in this package |
|---|---|
| Genus–environment associations | 12 fixed genera × 13 variables × 2 markers = 312 tests; full-composition CLR, locality-adjusted rank coefficients, 4,320 exhaustive permutations, BH156 and BY156 separately per marker |
| Alpha–environment associations | 2 indices × 13 variables × 2 markers = 52 tests; locality-adjusted rank coefficients, BH26 and BY26 separately per marker |
| Sensitivity | All 3,276 leave-one-unit-out coefficients and the specified pseudocount sensitivity; these do not add biological replicates |
| Alpha-category tests | Four fits independently recalculated over 180 distinct within-locality category allocations; current Holm families include the retained beta endpoints |
| Environmental category tests | All 28 rank-based partial-F comparisons recalculated; BH28, BY28 and Holm28; unblocked exact Kruskal–Wallis/BH28 retained as a separate sensitivity analysis |
| ANCOM-BC2 and dispersion | Current table identity, multiplicity and coefficient arithmetic verified; native models are not refitted by this numerical replay |
| Sequence processing, rarefaction and beta-distance fits | Outside this entry point; consult the sequence workflow and its source inputs |

The selected twelve genera are not renormalized as a separate composition. The CLR denominator uses all retained genus-rank lineage groups, including unresolved annotations. Both markers represent the same nine soil units; the one hundred ITS rarefaction draws are technical subsamples, not biological replication.

## What was actually verified

The current-only replay regenerated ten TSV files byte-for-byte identical to the corrected reference tables. The independent implementation passed 97 checks, including all 364 exact association permutation counts and the four alpha plus 28 environmental category fits. The maximum independent numerical discrepancy was 4.45 × 10^-16. Thirteen further checks cover retained ANCOM-BC2/dispersion arithmetic. The metadata preparation check passed 308 comparisons.

The executed reports are preserved in `validation/analysis/`. In a new replay, inspect `validation/REPLAY_COMPARISON.json`, `ASSOCIATION_INDEPENDENT_QA.json`, `RETAINED_STATISTICS_QA.json` and `METADATA_PREPARATION_QA.json`. `REFITTED_ALPHA_CATEGORY.tsv` and `REFITTED_ENVIRONMENT_CATEGORY.tsv` expose the recalculated category estimates. Their floating-point text formatting may differ from the original tables; the independent numerical check uses a documented tolerance of 1e-10.

`analysis/EXPECTED_TABLES.json` defines current association output hashes. Hash equality checks file identity, not the validity of field sampling, assays or biological interpretation.

## Interpretation retained

There are nine units, with six from Jaén (J) and three from San Ignacio (SI); low degradation occurs only in J. The models assume field independence and exchangeability within locality, which cannot be established from these files. Locality adjustment does not remove all potential vegetation/degradation confounding. No primary adjustment family yields P < 0.05; this is not evidence of equivalence or absence of an effect.

Metadata preparation records preserve the A1 particle-size normalization: the originally recorded values sum to 102%, and the analytical proportions were scaled to 100% and rounded to two decimals. This computational preparation is not a confirmed laboratory correction. Total elemental content was reported by the study author; units, digestion/assay details and detection-limit handling still require the laboratory source records.

The source tables and checks cannot establish contamination, toxicity, resistance, function or causality. Taxonomic genus assignments do not establish species identity or phenotype.
