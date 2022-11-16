```
site_name: GiniClust3
repo_url: https://github.com/rdong08/GiniClust3
site_description: Documents for GiniClust3
site_author: Rui Dong
pages:
    - Home: 'index.md'
    - Install: 'installation.md'
    - Tutorial: 'tutorials.md'
    - Functions:
        - 'calGini': 'func/calGini.md'
        - 'clusterGini': 'func/clusterGini.md'
        - 'calFano': 'func/calFano.md'
        - 'clusterFano': 'func/clusterFano.md'
        - 'generateMtilde': 'func/generateMtilde.md'
        - 'clusterMtilde': 'func/clusterMtilde.md'
        - 'plotGini': 'func/plotGini.md'
        - 'plotFano': 'func/plotFano.md'
    - About:
        - 'Release Notes': 'about/changelog.md'
        - 'License': 'about/license.md'
theme: 'readthedocs'

```
DISPbind is a DNA associated disorder protein analysis toolset.

Authors: Rui Dong (rdong@mgh.harvard.edu)

Maintainer: Rui Dong (rdong@mgh.harvard.edu)

## Features

* Genome mapping and bigwig generation of DisP-seq sequences ([Alignment](modules/align.md))
* Identification of DisP island ([Island](modules/island.md))

## Tutorial

We have developed a tutorial that demonstrates how DISPbind helps you to identify the genomeic binding of disorder protein and DisP island. Please check it below!

* [Installation and Setup](tutorial/setup.md)
* [Pipeline](tutorial/pipeline.md)
* [Mapping](tutorial/mapping.md)
* [Island](tutorial/island.md)

## Modules

DISPbind contains 3 modules. Each module functions as an independent component owning its distinctive duty. Meanwhile, they inteact with each other, and different circular RNA analysis pipelines are derived from different combinations of several modules. Understanding the detailed mechanism of each module could facilitate your circular RNA research.

List of Modules:

* [Align](modules/align.md)
* [Bam2bw](modules/bam2bw.md)
* [Island](modules/island.md)

## Citation

[Yu-Hang Xing\*, Dong R\*, Lukuo Lee, Shruthi Rengarajan, Nicolò Riggi, Gaylor Boulay and Miguel N. Rivera#.  *Under review*]

## License

Copyright (C) 2022 Rivera Lab. See [LICENSE](about/license.md) for license rights and limitations (MIT).
