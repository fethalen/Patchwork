Patchwork
---------

> **Note:** Patchwork is an unfinished product currently under development.

## Overview

Fragmented assemblies—generated from low-coverage sequence data or similar
technologies—can prove problematic in downstream analyses when used in a
phylogenomic context. Patchwork employs BLAST searches to "stitch" such
fragments into longer pieces. Once the orientation of the fragments has been
established, each fragment is also annotated, with the annotation possibly
being guided by a previously annotated genome.

## Installation

The following external programs must be installed:

* [DIAMOND](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/diamond/) (preferred) or [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)

Unless you already have a preferred method of installation, we recommend using
[Anaconda](https://www.anaconda.com/products/individual) in order to install
these programs. After Anaconda has been installed, the programs can be
installed by typing the following commands:

```
conda install -c bioconda mafft
conda install -c bioconda diamond
# Note that BLAST is not required if DIAMOND is installed
conda install -c bioconda blast
```

## Cite

Our manuscript is still in preparation, it will be posted here once a preprint
of the article is available.

© Animal Evolution and Biodiversity 2020
