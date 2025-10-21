FROM r-base:4.3.1

LABEL authors="Evan Biederstedt <evan.biederstedt@gmail.com>, Teng Gao <tgaoteng@gmail.com>" \
    version.image="0.1.3" \
    version.pagoda2="0.1.3" \
    description="r-base image to run numbat with R"

RUN apt-get update --yes && apt-get install --yes build-essential \
  libcurl4-gnutls-dev libxml2-dev libssl-dev libbz2-dev zlib1g-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libncurses5-dev libncursesw5-dev liblzma-dev libgit2-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libglpk-dev git autoconf gettext libtool automake  \
  samtools sudo

RUN apt-get update && apt-get install -y build-essential binutils gcc g++

RUN cd /usr/bin && \
  wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 && \
  tar -vxjf htslib-1.15.1.tar.bz2 && cd htslib-1.15.1 && make && sudo make install

RUN R -e 'chooseCRANmirror(ind=42); install.packages("BiocManager")'

RUN R -e 'chooseCRANmirror(ind=42); install.packages("ragg")'

RUN R -e 'chooseCRANmirror(ind=42); install.packages("pkgdown")'

RUN R -e 'chooseCRANmirror(ind=42); install.packages("devtools")'

RUN R -e 'devtools::install_github("YuLab-SMU/ggtree", dependencies=TRUE)'

RUN R -e 'devtools::install_github("kharchenkolab/numbat", dependencies=TRUE)'

RUN git clone https://github.com/kharchenkolab/numbat.git

RUN mkdir -p /tmp && chmod 777 /tmp

RUN chmod 777 /numbat/inst/bin/pileup_and_phase.R 


RUN mkdir data

# hg 38
RUN cd /data && wget -q https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz && gzip -d genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz 


## Note: not run; the user can follow these commands if they want hg19
# hg19
## RUN cd /data && wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz && gzip -d genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz


# hg38
RUN cd /data && wget -q http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip && unzip 1000G_hg38.zip 

# hg19
##RUN cd /data && wget http://pklab.med.harvard.edu/teng/data/1000G_hg19.zip && unzip 1000G_hg19.zip



RUN git clone https://github.com/single-cell-genetics/cellsnp-lite.git && cd cellsnp-lite && \
  autoreconf -iv && ./configure && make && sudo make install    

RUN wget -q https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz && cd .. && tar -xvzf Eagle_v2.4.1.tar.gz && cd /Eagle_v2.4.1 && cp eagle /usr/bin


WORKDIR /numbat


## NOTE: 
## Added 1000G Reference Panels for hg38  as follows:
## (1) Downloaded the file via wget in a Docker container
## # hg38
## wget http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip && unzip 1000G_hg38.zip
##
## (2) Then I save this container as an image using `docker container commit`
## cf https://docs.docker.com/engine/reference/commandline/push/
## i.e. `docker container commit c16378f943fe  pkharchenkolab/numbat-rbase:latest`
##
## Note: I didn't include hg19, though the procedure is the same using
## # hg19
## wget http://pklab.med.harvard.edu/teng/data/1000G_hg19.zip
## 

