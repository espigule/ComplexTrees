# Reproducibility

The article and poster compile using their included figures. A suitable TeX installation supplies the standard TeX fonts. The repository does not include font files.

Run `python3 tools/check_release.py` to verify all seven stored finite-exclusion benchmarks and run the finite exact-algebra tests. Their claims are explicitly bounded; these checks are not a new referee report or proof-assistant formalization.

To redraw selected stable surfaces, install the scientific Python requirements and run `python3 figure_sources/four_families/build_four_family_figures.py`. To rebuild the offline page, run `python3 interactive/build_interactive.py`. Positive interval cells and curve source data are retained; blank regions are not recoded as unstable.

The central atlas image comes from a finite DD survey. Reconstructing labels and specimen panels does not repeat or refine that parameter survey. The continuously developed full explorer remains at complextrees.com/4D/ and is separate from this frozen experiment.
