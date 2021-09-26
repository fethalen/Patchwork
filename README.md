> **Note:** Patchwork is an unfinished product, currently under construction 🚧.
<img src="https://github.com/fethalen/Patchwork/blob/main/patchwork_logo_500px.png" alt="Patchwork logo" width="225"/>

### Introduction

Fragmented genome assemblies—resulting from low-coverage sequencing or similar
technologies—can prove problematic in downstream analyses, when used in a
phylogenomic context. Patchwork employs local alignment searches to retrieve
homologous regions from a set of contigs and "stitch" them together. The
resulting homologs are directly suitable for use in phylogenomic studies.

### Graphical Overview

![Graphical Overview](https://github.com/fethalen/patchwork/blob/main/overview.png?raw=true)

### Features

* Align nucleotide sequences to one or more protein sequences
* Stitch overlapping or gappy sequences together based on a reference
* Find homologs, even in distantly-related taxa
* Written in [Julia](https://julialang.org/) and utilizing [DIAMOND](https://github.com/bbuchfink/diamond) for maximum speed 🐇

### Installation

> Note: We are working on a Conda build for Patchwork. In the future, 
> the user will be able to install the programming by running conda 
> install -c bioconda patchwork. Until then, please refer to 
> [these instructions](https://github.com/fethalen/Patchwork/wiki/4.-Installation)
> for installing Patchwork.

### Cite

Our manuscript is still in preparation, it will be posted here once a preprint
of the article is available.

© [Dept. for Animal Evolution and Biodiversity](https://www.uni-goettingen.de/en/80149.html) 2020
