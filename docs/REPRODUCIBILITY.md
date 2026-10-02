# Reproducibility

The manuscript and companion are in preparation. Their PDF and LaTeX sources are withheld pending author review. The submitted poster PDF is retained as the approved print artifact; the superseded poster sources have been removed from the current tree.

Run `python3 tools/check_release.py` to verify all seven stored finite-exclusion benchmarks and the finite exact-algebra checks. These are bounded scientific checks, not manuscript review or a proof-assistant formalization.

The existing `figure_sources/`, `computation/` and `verification/` records are unchanged. Install `requirements.txt` and run `python3 figure_sources/four_families/build_four_family_figures.py` to redraw selected stable surfaces. Run `python3 interactive/build_interactive.py` to rebuild the offline viewer. Positive interval cells retain their original scope; blank regions are not classified as unstable.

The full explorer at https://complextrees.com/4D/ develops separately from this preserved numerical experiment.
