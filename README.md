Patchwork
---------

> **Note:** Patchwork is an unfinished product, currently under construction 🚧.

### Introduction

Fragmented genome assemblies—resulting from low-coverage sequencing or similar
technologies—can prove problematic in downstream analyses, when used in a
phylogenomic context. Patchwork employs local alignment searches to retrieve
homologous regions from a set of contigs and "stitch" them together. The
resulting homologs are directly suitable for use in phylogenomic studies.

### Features

* Align nucleotide sequences to one or more protein sequences
* Stitch overlapping or gappy sequences together based on a reference
* Find homologs, even in distantly-related taxa
* Written in [Julia](https://julialang.org/) and utilizing [DIAMOND](https://github.com/bbuchfink/diamond) for maximum speed 🐇

### Graphical Overview

![Graphical Overview](https://github.com/fethalen/patchwork/blob/main/overview.png?raw=true)

### Installation

The sequence aligner [DIAMOND](https://github.com/bbuchfink/diamond) is required
to run Patchwork. We recommend using
[Anaconda](https://www.anaconda.com/products/individual) for installing this
program. After installing Anaconda, install DIAMOND by entering:

```bash
conda install -c bioconda diamond
```

### Cite

Our manuscript is still in preparation, it will be posted here once a preprint
of the article is available.

© [Dept. for Animal Evolution and Biodiversity](https://www.uni-goettingen.de/en/80149.html) 2020
