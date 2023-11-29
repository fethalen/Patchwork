<img src="https://github.com/fethalen/Patchwork/blob/main/patchwork_logo_500px.png" alt="Patchwork logo" width="225"/>

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)

Patchwork is an alignment-based program for retrieving and concatenating
phylogenetic markers from whole-genome sequencing (WGS) data. The program
searches the provided DNA query contigs against one or more amino acid reference
sequences. Multiple, overlapping hits are merged to derive a single, continuous
sequence for each reference sequence.

### Features

* Align nucleotide sequences to one or more protein sequences
* Works with already assembled contigs _or_ raw reads
* Stitch overlapping or gappy sequences together based on a reference
* Find homologs, even in distantly-related taxa
* 🐇 Written in [Julia](https://julialang.org/) and utilizing
  [DIAMOND](https://github.com/bbuchfink/diamond) for maximum speed

### Graphical Overview

![Graphical Overview](https://github.com/fethalen/patchwork/blob/main/overview.png?raw=true)

### Quick installation

We are currently working on a Conda build. In the future,
the user will be able to install this program by running `conda
install -c bioconda patchwork`. Until then, please refer to
[these instructions](https://github.com/fethalen/Patchwork/wiki/4.-Installation)
for installing from source. It is now also possible to [install
Patchwork using Docker](https://github.com/fethalen/Patchwork/wiki/4.-Installation#installing-patchwork-with-docker).

### Documentation

Please see our [Wiki](https://github.com/fethalen/Patchwork/wiki).

### Cite

Thalen, F., Koehne, C. G., & Bleidorn, C. (2022). Patchwork: alignment-based retrieval and concatenation of phylogenetic markers from genomic data. bioRxiv.

BibTeX-format:

```bibtex
@article{thalen2022patchwork,
  title={Patchwork: alignment-based retrieval and concatenation of phylogenetic markers from genomic data},
  author={Thalen, Felix and Koehne, Clara Gisela and Bleidorn, Christoph},
  journal={bioRxiv},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}
```

© [Dept. for Animal Evolution and Biodiversity](https://www.uni-goettingen.de/en/80149.html) 2020
