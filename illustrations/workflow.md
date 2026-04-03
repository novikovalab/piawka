---
layout: default
title: Illustrations guide
parent: Technical details
nav_exclude: true
---

# Illustration workflow guide

This page describes how to produce hand-drawn-style illustrations for the piawka documentation website, consistent with the creature-based logo style.

---

## Target figures

The following figures are proposed for the documentation. Each is described by location, content, and recommended dimensions.

### 1. Overview page — workflow diagram

**Location**: `index.md`, after the subcommands table  
**Content**: Shows the overall `piawka` pipeline:

```
VCF (.vcf.gz) ──→ piawka calc ──→ BED output
                         │
                         ├──→ piawka sum  ──→ summarized BED
                         ├──→ piawka filt ──→ filtered BED
                         └──→ piawka dist ──→ distance matrix
```

**Style**: leech mascot dragging a VCF file on the left, output documents on the right, connected by arrows.

---

### 2. Subcommands page — calc illustration

**Location**: `subcommands.md`, calc section  
**Content**: Shows a VCF being queried by `tabix` per region, with allele counts pooled per group.

---

### 3. Statistics page — pi illustration

**Location**: `statistics.md`, pi section  
**Content**: Two alleles drawn from a pool; probability of difference indicated by an arrow.

---

### 4. Technical page — SNP retrieval

**Location**: `technical.md`, SNP retrieval section  
**Content**: Three population bubbles with allele counts. Full cohort = multiallelic; per-group = biallelic.

---

### 5. Technical page — increment/finalize framework

**Location**: `technical.md`, statistics framework section  
**Content**: Flowchart: VCF site → increment → accumulate num/den → finalize → print.

---

### 6. Technical page — parallelism

**Location**: `technical.md`, parallel processing section  
**Content**: Parent + N children with FIFO buffers and tmp files.

---

## Recommended tools

### Option A: Excalidraw (recommended)

[Excalidraw](https://excalidraw.com) is a free, browser-based tool that produces hand-drawn-style SVG diagrams natively. It supports:

- Export as `.svg` (scalable, embeddable)
- Creature/sketch style that matches the piawka logo aesthetic
- Collaborative editing via share links

**Workflow:**

1. Open [excalidraw.com](https://excalidraw.com).
2. Draw the diagram using the freehand and shape tools.
3. Export as SVG: `⋮ → Export image → SVG`.
4. Save to `assets/img/<figure-name>.svg` in the `gh-pages` branch.
5. Embed in Markdown: `![Description](assets/img/<figure-name>.svg)`.

---

### Option B: Inkscape

[Inkscape](https://inkscape.org) gives full SVG control and supports the "sketch" filter that mimics hand-drawn lines.

**To apply a hand-drawn look:**

1. Draw shapes in Inkscape.
2. Select all → `Filters → Distort → Rough and Jagged` (or similar).
3. Export as **Plain SVG** for small file size.

---

### Option C: Python (matplotlib + Rough.js)

For data-driven diagrams (e.g., allele frequency spectra), generate figures with `matplotlib` and then apply a hand-drawn style with the `mplfinance` rough patch or save as SVG and post-process with [rough.js](https://roughjs.com).

---

## SVG embedding

Once SVG files are placed in `assets/img/`, embed them in any page:

```markdown
![Workflow diagram]({{ site.baseurl }}/assets/img/workflow.svg)
```

For figures with captions, use HTML directly:

```html
<figure>
  <img src="{{ site.baseurl }}/assets/img/workflow.svg"
       alt="piawka workflow" style="max-width:600px">
  <figcaption>Fig. 1 — Overall piawka pipeline.</figcaption>
</figure>
```

---

## Style guidelines

- **Line weight**: 1.5–2px strokes, slightly irregular (hand-drawn feel)
- **Colours**: use the piawka palette (`#FA9E89` for accents, `#3D3D3D` for text/lines, white fill)
- **Font**: a rounded, hand-written-style font (e.g., "Caveat" from Google Fonts)
- **Size**: keep SVG viewBox at 600–800px wide; the site layout constrains max-width anyway
- **Mascot**: use the leech character from `logo/logo.svg` as a reference; resize and recolour freely
