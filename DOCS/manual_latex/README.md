# GHATS LaTeX Manual Source

This directory contains the LaTeX version of the GHATS manual.

The intended workflow is:

1. Keep the legacy Pages, Word, and PDF manual files in `DOCS/legacy_manual_sources/`.
2. Use the Word document as the main text source when recovering old text.
3. Use `GHATS_Manual_v3.3.0.pdf` as the visual/style reference.
4. Review and update the manual one chapter at a time against the current routines.

The current LaTeX manual has chapter and appendix sources, a reconstructed
title page, selected walkthroughs, and a generated `main.pdf`.

Build the distributed PDF with XeLaTeX:

```sh
xelatex main.tex
xelatex main.tex
```

XeLaTeX is preferred because it uses the thin sans-serif title-page font
configured in `main.tex`, which is closer to the original Pages/PDF manual.
