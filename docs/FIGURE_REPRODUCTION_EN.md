# Reproducing the current scientific figures

This component contains the publication figures and rendering inputs for the corrected Yanayacu edition. Files in `figures/MAIN` and `figures/SUPPLEMENTARY` retain the source PDF, SVG and PNG bytes.

## Run in a fresh output directory

From the public release root:

```bash
python -m pip install -r figure_code/requirements-render.txt
python figure_code/reproduce_figures.py --output ../yanayacu_figures_replay
```

The output directory must not already exist. The command does not overwrite the distributed figures or refit statistical tests. It checks source hashes, renders into the new directory, writes a gallery and compares the saved PDF numbers and dimensions against the distributed current figures. Numerical input matrices remain unmodified.

## Exact scope

| Operation | Figures | Meaning |
|---|---|---|
| Rebuild from frozen numerical tables | F02, F03, F05, F06, S03, S04, S05, S06, S07, S16, S17, S21, S23 | Thirteen plots are rendered from sample measurements, coordinates and test results. S07 recalculates its descriptive correlation display; no inferential model is refitted. |
| Native vector artwork plus source-based text update | S09, S11 | Original quantitative geometry is preserved. S09 receives the current test probabilities; S11 receives the clarified technical-sensitivity title. |
| Preserve the current publication artifacts byte-for-byte | F01, F04, S01, S02, S08, S10, S12, S13, S14, S15, S18, S19, S20, S22, S24, S25 | These sixteen figure sets remain available as PDF/SVG/PNG and are copied by the figure replay. The replay does not claim to rebuild their original statistical or sequence-processing calculations. |

`data/figure_inputs/native/` contains the two vector geometry inputs used by the builders. Numeric tables and coordinates are in `data/figure_inputs/`. The hashes and original relative source locations are in `validation/figures/SOURCE_INVENTORY.tsv`.

## Source-to-figure mapping

* `build_associations.py`: four current association/sample tables feed the eleven association and soil-correlation figures.
* `build_alpha.py`: the alpha measurements, technical-subsampling summaries and current category-test table feed F02.
* `build_beta.py`: frozen PCoA coordinates, axis fractions, dispersion distances and current category probabilities feed F03. The current four-endpoint ITS category family is read directly; the older two-endpoint category annotation is not used. Dispersion remains a separate diagnostic family.
* `annotate_s09.py`: reads the current category-association table to update the native S09 footer.
* `clarify_s11_title.py`: updates the technical title while retaining all quantitative artwork.
* `FIGURE_CAPTIONS.json`: authoritative current English legends for every figure.

## Executed verification

The standalone replay was executed on 26 November 2026 with Python 3.13.12, NumPy 2.4.6, pandas 3.0.6, Matplotlib 3.10.9, Pillow 12.2.0 and PyMuPDF 1.27.2.2. These versions describe the rendering check; they do not identify the original model-fitting environments. The original current-edition rendering environment is retained in `data/figure_inputs/SOURCE_RENDER_ENVIRONMENT.json`.

All 31 PDF numerical-token sequences and page dimensions matched the current-source figures. All 2,272 association plot objects, all 36 alpha-source rows and all 36 beta-source rows matched their current numerical sources. The 93 distributed publication artifacts are byte-identical to the current edition. Sixteen replayed figure sets are byte-identical; rebuilt plots can differ at the pixel or file-byte level because of library/font versions. These differences are measured explicitly in `validation/figures/FIGURE_REPRODUCTION_CHECK.json`. The rebuilt figures are therefore assessed separately from the byte-identical distributed assets.

The checks concern source identity, graphical quantities, numerical labels and rendering. They do not independently establish sampling independence, laboratory assay validity, statistical model adequacy or scientific significance. See the analysis reproduction guide for numerical-model scope.
