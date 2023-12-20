FROM ubuntu:20.04 
ENV http_proxy=http://proxy.nhri.org.tw:3128 
ENV https_proxy=http://proxy.nhri.org.tw:3128 
ENV TZ=Asia/Taipei 
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone 
RUN apt-get update && apt-get install -y \ 
    python3 \ 
    python3-pip \ 
    python3-setuptools \ 
    python3-dev \ 
    libbz2-dev \ 
    liblzma-dev \ 
    libncurses5-dev \ 
    libncursesw5-dev \ 
    zlib1g-dev \ 
    libcurl4-gnutls-dev libssl-dev \ 
    cmake unzip git wget libz-dev vim autoconf curl 
RUN pip3 install biopython pandas lxml six ete3 openpyxl

#Download nanoBSI
WORKDIR /opt
RUN git clone https://github.com/jade-nhri/nanoBSI.git
WORKDIR /opt/nanoBSI
RUN chmod +x *.py

#Install Minimap2 2.20-r1061 
WORKDIR /opt 
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.20/minimap2-2.20_x64-linux.tar.bz2 | tar -jxvf - 
RUN mv minimap2-2.*/ minimap2 

#Install blast2.15.0
WORKDIR /opt
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN tar zxvf ncbi-blast-2.15.0+-x64-linux.tar.gz


#Install ResFinder 4.4.2
WORKDIR /opt 
RUN pip3 install resfinder
RUN git clone https://bitbucket.org/genomicepidemiology/resfinder_db/
RUN git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/
RUN git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/
RUN git clone https://bitbucket.org/genomicepidemiology/kma.git
ENV CGE_RESFINDER_RESGENE_DB="/opt/resfinder_db"
ENV CGE_RESFINDER_RESPOINT_DB="/opt/pointfinder_db"
ENV CGE_DISINFINDER_DB="/opt/disinfinder_db"
WORKDIR /opt/kma
RUN make
WORKDIR /opt/pointfinder_db
RUN python3 INSTALL.py /opt/kma/kma_index

#Install Centrifuge 1.0.4
WORKDIR /opt
RUN git clone https://github.com/infphilo/centrifuge
WORKDIR /opt/centrifuge
RUN make
RUN make install prefix=/usr/local

#Install Flye, 2.9.3
WORKDIR /opt
RUN git clone https://github.com/fenderglass/Flye.git
WORKDIR /opt/Flye
RUN make

#Install seqkit 2.6.1
WORKDIR /opt 
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.6.1/seqkit_linux_amd64.tar.gz
RUN tar -zxvf seqkit_linux_amd64.tar.gz 
RUN cp seqkit /usr/local/bin 


WORKDIR / 
#set path 
ENV PATH $PATH:/opt:/opt/:/opt/minimap2/:/opt/centrifuge:/opt/Flye/bin/:/opt/ncbi-blast-2.11.0+/bin:/opt/nanoBSI/

