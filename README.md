# nanoBSI
Comprehensive pathogen identification and antimicrobial resistance prediction from positive blood cultures using nanopore sequencing technology

**To run with Docker**

``git clone https://github.com/jade-nhri/nanoBSI.git``

``cd nanoBSI``

``docker build -t "nanobsi:v1" ./``

``docker run -h nanobsi --name nanobsi -i -t -v /:/MyData nanobsi:v1 /bin/bash``


Installation
------------
**Installation from source**

``cd /opt``

``git clone https://github.com/jade-nhri/nanoBSI.git``

``cd nanoBSI``

``chmod +x *.py``

``export PATH="$PATH:/opt/nanoBSI"``


## Dependencies


- [minimap2-2.20](https://github.com/lh3/minimap2)
- [seqkit-v2.6.1](https://github.com/shenwei356/seqkit)
- [blast-2.15.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
- [Flye-2.6.1](https://github.com/fenderglass/Flye)
- [Centrifuge-1.0.4](https://github.com/infphilo/centrifuge)
- [ResFinder-4.4.2](https://pypi.org/project/resfinder/)

## FASTQ files

https://www.dropbox.com/scl/fi/lbdrmes0kl8954xdmhwe6/fastq_pass.tar.gz?rlkey=y43fxw50mbcedk4plappv7y8d&dl=0

## Quick start

- [Readme](https://www.dropbox.com/scl/fi/mtpmd84ax8zh216hz4a6w/Manual-of-nanoBSI.pdf?rlkey=8p2oabkj7b2px4rg95uj6fa9s&dl=0)
