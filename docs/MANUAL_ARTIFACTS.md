# Manually produced figure artifacts

Some of what appears in the manuscript was not produced by this pipeline. This
file records exactly which parts, so nobody spends time trying to reproduce a
hand-adjusted layout from code.

Source: `archive/final_package_v3_original/04_METADATA/manual_editing_provenance.tsv`.

## Affinity — multi-panel assembly

All multi-panel figures were composed and cosmetically edited in Affinity.

**What the pipeline guarantees:** the numerical content, the selected features,
and each source panel as its own PDF/PNG pair.

**What it does not:** the page composition, panel lettering, cropping, shared
legends, or any typographic adjustment.

The audit's own wording, which is the right standard:

> the analytical scripts should be expected to reproduce the underlying source
> plots/data, not necessarily the exact final page composition.

**To close this gap:** deposit the Affinity sources (`.afdesign`) or exported
SVG/PDF of the final compositions.

## Cytoscape v3.10.3 — network layouts

Figures 6B-F and S5A-B are Cytoscape renderings whose layout was adjusted by
hand.

**What the pipeline guarantees:** network membership, every node and edge
attribute, and import-safe tables:

- `R24_ceRNA_cytoscape_export` → `ALL_MFE20_nodes.csv`, `ALL_MFE20_edges.csv`,
  per-module node/edge/triplet tables, a visual-style template, and an export manifest
- `R25_cytoscape_key_validation` → `*_CYTOSCAPE_SAFE.tsv`, `ALL_MFE20_edges_MINIMAL.tsv`,
  and an ID-match report

**What it does not:** node coordinates. Re-importing these tables reproduces the
same graph with a different layout.

**To close this gap:** deposit the saved `.cys` session, or export the final node
x/y coordinates plus the visual style.

## Figure 1 and Figure 2A

| Panel | Content | Source |
|---|---|---|
| Figure 1 A-D | Experimental workflow, immunofluorescence, area/perimeter | Experimental data, manual assembly |
| Figure 2 A | Integrated bioinformatic workflow | Hand-drawn diagram |
| Figure S1 A-E | Phenotypic, qPCR and western blot validation | Experimental data |

None of these is transcriptomic output. Preserve the editable diagram source for
Figure 2A if it exists.

## A note on the visualization notebook

`Babi analysis/nrvm_transcriptome_visualization.qmd` contains more than forty
characterization panels and **zero `ggsave()` calls** — every one of them exists
only inside the rendered HTML. If any of those panels is wanted for the
supplement, it has to be re-emitted through `save_figure()` so there is a file
on disk. Per the decision taken for this rebuild, that notebook is not in the
repository; the originals are in `SUZ2026/Babi analysis/`.
