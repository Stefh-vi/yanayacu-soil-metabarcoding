# Table legends

## Table 1. Sequence accounting through the 16S and ITS workflows

Counts are pooled across nine soil units per marker. The same soil units contribute both markers. Units change from paired reads to forward reads or reconstructed sequences at the indicated stages; these counts must not be interpreted interchangeably. Taxonomically retained sequences comprise the prokaryotic 16S and explicitly fungal ITS libraries. Rarefaction totals refer to one draw across all nine units. The final percentage is a sequence-to-raw-pair count ratio, not a fraction of individual reads. ITS fixed-length truncation was not applied.

## Table 2. Diversity associations with degradation category after locality adjustment

Category tests use the same nine soil units, with three per category, and condition on locality. Effect and residual degrees of freedom are 2 and 5. Raw P values enumerate 4,320 within-locality arrangements; Holm adjustment covers three endpoints for 16S and four for ITS. Total R² and partial R² after locality adjustment have different denominators. Alpha estimates are marker-specific; the two ITS distances are alternative descriptions of the same community. PERMDISP results are reported separately in Supplementary Table S06.

## Supplementary Table S01. Sequence counts for individual soil libraries

Counts are provided for each sample, marker and processing stage, preserving the original units: paired reads, forward reads or reconstructed sequences. Stage order defines the workflow sequence. Each marker represents nine soil units rather than an independent set of biological replicates. Counts at successive stages with different units cannot be subtracted as read losses. The table provides the sample-level values underlying the totals in Table 1; no statistical test is associated with these accounting data.

## Supplementary Table S02. Changes in sequence counts between processing stages

Counts, count ratios relative to raw read pairs and losses between comparable consecutive stages are listed for each sample and marker. Percent decreases use the previous-stage count only when the counting units match. The comparison-type and interpretation fields identify changes from read pairs to forward reads or reconstructed sequences. Not applicable (n.a.) denotes the initial input or an invalid subtraction between different counting units, rather than a missing biological observation. Ratios involving distinct units are explicitly descriptive accounting ratios.

## Supplementary Table S03. Primary alpha-diversity estimates

Observed ASV richness and Shannon diversity, using natural logarithms, are reported for each of the nine soil units and each marker. The 16S estimates use one rarefaction draw at 174,661 retained prokaryotic sequences per unit; ITS estimates are means of 100 draws at 450,343 fungal sequences per unit. Locality and degradation category accompany each observation. These are the primary estimates used in Figure 2 and category-association tests. Technical subsampling summaries are reported separately in Supplementary Table S04.

## Supplementary Table S04. Technical variation in alpha diversity under rarefaction

Means and 2.5th–97.5th percentiles of observed ASV richness and Shannon diversity from 100 subsamples without replacement at the marker-specific depth for every soil unit. Percentile ranges describe technical subsampling variation, not confidence intervals for biological variation. The relation-to-primary-estimate field distinguishes the primary single-draw 16S estimate from its additional technical summary; ITS uses the 100-draw mean as its primary estimate. The number of biological units remains nine for each marker. No significance tests are attached to these percentile summaries.

## Supplementary Table S05. Alpha-diversity associations with environmental variables

Post hoc associations of two alpha indices with 13 variables, totaling 26 tests per marker. Coefficients correlate globally ranked alpha and environmental values after removing locality means across nine soil units. Two-sided raw P values use all 4,320 within-locality arrangements, including the observed configuration. Benjamini–Hochberg (BH) adjustment is primary; Benjamini–Yekutieli (BY) adjustment is a sensitivity analysis, each across 26 tests per marker. Fields prefixed q contain these adjusted P values, not Storey q values. Leave-one-out coefficient ranges describe sensitivity and are not confidence intervals.

## Supplementary Table S06. Multivariate dispersion diagnostics

PERMDISP compares distances to category spatial medians, with small-sample bias adjustment, among three degradation categories in nine soil units. Effect and residual degrees of freedom are 2 and 6. Permutations are restricted within locality, but the dispersion models are not conditioned on locality. The 16S Bray–Curtis diagnostic has an unadjusted P value; the two ITS diagnostics use a separate Holm family of size two. These tests address dispersion rather than category centroids. Not-applicable R² fields do not denote failed category tests.

## Supplementary Table S07. Distances used in Bray–Curtis dispersion diagnostics

Distance from each soil unit to its degradation-category spatial median, with the recorded small-sample bias adjustment. Values correspond to Bray–Curtis dissimilarities from rarefied 16S counts or full ITS fungal-library proportions. The table contains nine observations per marker and supplies Figure 3c–d. Category and locality identify the sampling structure. These distances are inputs to the PERMDISP diagnostics in Supplementary Table S06; they are neither confidence limits nor independent replicates of their originating soil units.

## Supplementary Table S08. NMDS coordinates and pointwise stress contributions

Two-dimensional coordinates and pointwise stress contributions for all nine soil units in each marker. Category and locality labels accompany the coordinates. NMDS positions are descriptive representations of the primary Bray–Curtis matrices and have no explained-variance percentages. Nearly coincident coordinates are retained without jitter. Coordinates and stress contributions do not increase the biological sample size. The fit diagnostics in Supplementary Table S09 are needed to interpret the near-zero-stress, nearly one-dimensional configurations shown in Supplementary Figure S14.

## Supplementary Table S09. NMDS fitting diagnostics

One record per marker documents sample size, dimensionality, random seed, stress-1, random starts, repeated best solutions, engine iterations and stopping information. Fit warnings and singular-value ratios characterize the nearly one-dimensional configurations. Both analyses used nine soil units. Near-zero stress and repeated best solutions do not establish that a two-dimensional ordination provides an informative ecological representation. Nonmetric and linear fit R² describe distance-fit diagnostics rather than variation explained by degradation category. The primary community representation remains PCoA.

## Supplementary Table S10. Pairwise distances for NMDS Shepard diagnostics

Observed Bray–Curtis dissimilarities, native ordination distances and monotonic fitted distances for all 36 distinct sample pairs per marker. The 72 rows represent two descriptions of the same nine biological units, not 72 independent observations. Distances preserve the native solver scale used for the Shepard plots in Supplementary Figure S15. Uniform rescaling of site coordinates for display does not change stress. No smoothing, pair exclusion or additional hypothesis test was applied to these diagnostic data.

## Supplementary Table S11. Genus-level annotation coverage

Counts and fractions assigned to each annotation category in every retained library. Fractions range from zero to one and use the complete retained prokaryotic or fungal-library sequence count as denominator. The definition field preserves the marker-specific distinction between checked names, other candidate names and unresolved assignments for 16S, versus nominal genus-path annotations and unresolved assignments for ITS. These criteria do not establish equivalent certainty between markers or validate sequence identity. All categories, including unresolved assignments, contribute to Supplementary Figure S01.

## Supplementary Table S12. Associations of selected genera with environmental variables

Post hoc associations of 12 selected genus-rank lineages with 13 variables, totaling 156 tests per marker. Coefficients correlate globally ranked centered log-ratio abundances and environmental values after removing locality means across nine units. The primary CLR pseudocount is 0.5; all 1,143 16S or 1,185 ITS lineage groups form the reference. Two-sided raw P values use 4,320 within-locality arrangements. BH and BY adjustments each cover 156 tests per marker; q-prefixed fields denote adjusted P values. Leave-one-out ranges and the one-count pseudocount comparison quantify sensitivity rather than confidence intervals.

## Supplementary Table S13. Genus abundances and corresponding environmental observations

Sample-level counts, complete retained-library denominators, relative abundances, centered log-ratio (CLR) values and 13 environmental measurements for the 12 selected lineages per marker. Relative-abundance fractions and percentages are provided separately and are not renormalized across selected taxa. CLR values use pseudocounts of 0.5 and one count for primary and sensitivity calculations. Repeated sample rows across lineages do not create additional biological replicates: there are nine shared soil units. These observations underlie the coefficients in Supplementary Table S12.

## Supplementary Table S14. Covariation among environmental variables

All 78 distinct pairs among the 13 environmental variables, using the same nine soil units. Raw Spearman coefficients are reported alongside correlations of global midranks after removal of locality means. These descriptive coefficients do not carry P values or a multiplicity family. They summarize shared variation relevant to interpreting separate genus–environment and alpha–environment associations. High correlation does not establish that quantities were measured independently or derived from one another. Rank coefficients do not resolve missing laboratory-unit or analytical-method documentation.

## Supplementary Table S15. ASV intersections among degradation categories

Mutually exclusive intersections of category-specific ASV sets for each marker. An ASV qualifies for a category when detected in at least two of its three soil units in the complete retained, unrarefied table. Columns identify the categories in each intersection, counts of intersecting ASVs and the detection rule. ASVs excluded from all category sets are reported separately. These depth-sensitive inventories are descriptive rather than significance tests. The 16S and ITS ASV sets are separate marker-specific entities and are not intersected with each other.

## Supplementary Table S16. Metadata for the nine soil units

Sample identifiers, degradation category, locality, vegetation, texture and the recorded environmental measurements used in the analyses. J denotes Jaén and SI San Ignacio; low degradation has no San Ignacio unit. Texture fractions are percentages and pH is dimensionless. Other units are explicitly unconfirmed because primary laboratory documentation was unavailable. Numeric values are retained on their recorded scales rather than assigned assumed units. The fertile-soil sample is outside this nine-unit analytical cohort. Data transformations and source provenance are documented separately in the reproducibility records.

## Supplementary Table S17. ANCOM-BC2 global tests of degradation associations

Global tests from the model locality + degradation category for all eligible genus-rank lineage groups: 356 for 16S and 179 for ITS. Holm adjustment is performed within each marker, including unresolved lineages. Abundance, detection, design-rank, convergence, endpoint-validity and pseudocount-sensitivity fields accompany the test statistics. Auxiliary coefficient failures must be distinguished from valid primary degradation endpoints. “None reported” in a failure field is a textual encoding, not proof of biological validity. These post hoc model-based tests use nine soil units and do not directly measure absolute abundance.

## Supplementary Table S18. ANCOM-BC2 degradation contrast estimates

Natural-log coefficients for moderate and high degradation relative to low degradation, adjusted for locality, across the same nine soil units. Estimates, standard errors and coefficient ± one standard error are reported; these ranges are not 95% confidence intervals. The global raw and Holm-adjusted P values are repeated for reference and are not separate pairwise P values. Global Holm families contain 356 16S and 179 ITS lineages. Validity, convergence, auxiliary-failure and sensitivity fields retain the interpretation limits of each fit. The table supplies Supplementary Figures S18–S19.

## Supplementary Table S19. Phylum-level percentages in 16S libraries

Sequence percentages for each of nine retained prokaryotic libraries, with the complete retained count as denominator. The ten leading resolved phylum labels are displayed individually, remaining labels form Other, and unresolved phylum is retained. Percentages include all three components and are not renormalized after selecting leading labels. This table contains the numerical values underlying Figure 4a, including segments too small for printed percentage labels. Sequence proportions describe relative taxonomic composition rather than absolute microbial abundance; no category significance test is included.

## Supplementary Table S20. Phylum-level percentages in ITS libraries

Sequence percentages for each of nine explicitly fungal libraries, with the complete retained fungal count as denominator. The ten leading resolved phylum labels are displayed individually, remaining labels form Other, and unresolved phylum is retained. Percentages include all three components and are not renormalized after selecting leading labels. This table contains the numerical values underlying Figure 4b, including segments too small for printed percentage labels. Sequence proportions describe relative taxonomic composition rather than absolute fungal abundance; no category significance test is included.

## Supplementary Table S21. Environmental associations with degradation category

Tests for all 28 recorded numeric environmental variables across nine soil units. Primary partial F tests compare rank(value) ~ locality + degradation with the locality-only model, using effect/residual degrees of freedom of 2/5. Raw conditional P values enumerate 4,320 within-locality arrangements, representing 180 distinct category allocations. BH adjustment across 28 variables is primary; BY and Holm adjustments are sensitivity summaries. Partial R² refers to ranks after locality adjustment. Exact unblocked Kruskal–Wallis tests across 1,680 allocations, with their own BH adjustment, are separate sensitivity analyses and do not replace the primary model.

## Supplementary Table S22. Environmental summaries by degradation category

Number of observations, median, minimum, maximum, arithmetic mean and sample standard deviation for each of 28 environmental variables within each degradation category. Each group contains three soil units; the table therefore has 84 variable-by-category rows. Values preserve the recorded scale. Standard deviations describe the observed units and are not standard errors or confidence intervals. Except for dimensionless pH and texture percentages, measurement units require laboratory confirmation. Category hypothesis tests and locality adjustment are provided in Supplementary Table S21 rather than inferred from these descriptive summaries.

## Supplementary Table S23. Differences between category medians for environmental variables

Observed category medians and their pairwise differences for each of 28 environmental variables. The reported difference is the second-category median minus the first-category median, with the direction specified in separate columns. Each category contains three soil units. These 84 contrasts are descriptive and neither post hoc significance tests nor confidence intervals. No significance letters or pairwise P values are implied. Units follow the recorded scales and their documented uncertainty. Locality-adjusted category tests are reported in Supplementary Table S21.

## Supplementary Table S24. Selection of genus-rank lineages for environmental associations

Twelve lineages per marker, with full ancestral-path identifiers, display order, prevalence, mean relative abundance, member-ASV counts and selection criteria. Eligibility required detection in at least six of nine units and mean retained-library abundance of at least 0.1%; displayed lineages were selected by abundance rather than association P values. Full-library fractions and percentages are both provided. The selected taxa are shared across the abundance and environmental-association figures. Nominal genus labels are reference annotations rather than validated species identities.

## Supplementary Table S25. Complete retained-library denominators

Retained sequence counts for each of the nine soil units in each marker, together with the size of the complete genus-rank reference: 1,143 lineages for 16S and 1,185 for ITS. Counts define the denominators for reported sequence proportions, including selected-genera displays. They are not restricted to the 12 displayed lineages and are not rarefied denominators. The complete reference also defines the centered log-ratio transformation. The 18 rows describe two markers from nine shared biological units.

## Supplementary Table S26. Alpha diversity and corresponding environmental observations

Primary observed-ASV and Shannon estimates, sampling metadata and all 13 environmental values for both markers. The table records the estimator, number of technical draws and rarefaction depth for every observation. Its 36 rows represent two indices × two markers × nine shared soil units, not 36 biological replicates. The 16S estimator is one draw at 174,661 sequences; the ITS estimator is the mean of 100 draws at 450,343 sequences. These observations underlie the associations in Supplementary Table S05.

## Supplementary Table S27. Rank dependence among nitrogen, organic matter and carbon

Raw Spearman correlations for the three distinct pairs among recorded total nitrogen, organic matter and total carbon in nine soil units. These coefficients are descriptive and have no attached P values or multiplicity family. Strong dependence means that concordant associations cannot be treated as three independent lines of environmental evidence. Correlation alone does not show whether the laboratory measured these quantities independently or calculated one from another. The measurements remain in the analyses; their provenance and laboratory methods require primary documentation.
