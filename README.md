<img src="https://github.com/fethalen/Patchwork/blob/main/patchwork_logo_500px.png" alt="Patchwork logo" width="225"/>

Patchwork is an alignment-based program for retrieving and concatenating phylogenetic markers from whole-genome sequencing (WGS) data. The program searches the provided DNA query contigs against one or more amino acid reference sequences. Multiple, overlapping hits are merged to derive a single, continuous sequence for each reference sequence.

### Features

* Align nucleotide sequences to one or more protein sequences
* Stitch overlapping or gappy sequences together based on a reference
* Find homologs, even in distantly-related taxa
* Written in [Julia](https://julialang.org/) and utilizing [DIAMOND](https://github.com/bbuchfink/diamond) for maximum speed 🐇

### Graphical Overview

![Graphical Overview](https://github.com/fethalen/patchwork/blob/main/overview.png?raw=true)

### Quick installation

We are currently working on a Conda build. In the future, 
the user will be able to install the programming by running `conda 
install -c bioconda patchwork`. Until then, please refer to 
[these instructions](https://github.com/fethalen/Patchwork/wiki/4.-Installation)
for installing from source. It is now also possible to [install
Patchwork using Docker](https://github.com/fethalen/Patchwork/wiki/4.-Installation#installing-patchwork-with-docker).

### Documentation

Please see our [Wiki](https://github.com/fethalen/Patchwork/wiki). 

### Cite

Our manuscript is still in preparation, it will be posted here once a preprint
of the article is available.

© [Dept. for Animal Evolution and Biodiversity](https://www.uni-goettingen.de/en/80149.html) 2020
