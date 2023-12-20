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
