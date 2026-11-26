# Native software and source notices

Native software output, copied package function definitions and project execution scripts have distinct provenance. No blanket project license overrides third-party software licenses.

- **ANCOMBC 2.14.0**: the retained `ANCOMBC_DESCRIPTION.txt` records **Artistic-2.0**. `SOURCE_NATIVE_*.txt` and diagnostic/instrumented function snapshots derive from that installed package. Copyright, authorship and license remain with the package authors. The execution code records its diagnostic instrumentation and preserves the original function snapshots. See the [Bioconductor ANCOMBC package](https://bioconductor.org/packages/ANCOMBC/) and [Artistic License 2.0](https://opensource.org/license/artistic-2-0).
- **vegan 2.6.8 / permute**: native model objects, installed-version records and project R calls are retained as statistical provenance. The project scripts call the installed software; they do not replace its license. See the [vegan source repository](https://github.com/vegandevs/vegan) and its package license information.
- **MultiQC / FastQC / Cutadapt / QIIME 2 / DADA2**: reports, command records and environment specifications retain their software attribution. Their respective licenses govern third-party components, embedded report assets and source snippets. The MultiQC HTML notices remain in the included report copies.

The numerical-extension ANCOM source modifies only in-memory progress/diagnostic instrumentation, as documented in the execution source. The original mathematical function snapshots and final native objects are included. No newly authored package-level license is assigned here to copied third-party code.
