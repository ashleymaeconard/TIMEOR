FROM rocker/r-ver:3.6.1

RUN apt-get update && apt-get install -y \
    libjpeg-dev \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    libxml2-dev \
    libssl-dev \
    wget \
    python-pip \
    python2.7 

# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    sudo R -e 'install.packages(c("openssl", "httr", "rvest", "xml2"), dependencies=TRUE)' && \
    sudo R -e 'install.packages(c("UpSetR, corrplot","autoplotly","BiocManager","broom","cluster","d3heatmap","data.table","devtools","doMC","doRNG","dplyr","DT","europepmc","factoextra","farver","fastmatch","fpc","ggforce","ggfortify","ggplot2","ggplotify","ggridges","GlobalOptions","graphlayouts","gridGraphics","heatmaply","markdown","plotly","png","polyclip","promises","reticulate","rmarkdown","rvcheck","shiny","shinyalert","shinycssloaders","shinydashboard","shinydashboardPlus","shinyjs","shinyLP","shinyWidgets","stringr","tibble","tidygraph","tidyr","tidyverse","triebeard","tweenr","urltools","vegan","zoo"))' && \
    sudo R -e 'BiocManager::install("clusterProfiler")' && \
    sudo R -e 'BiocManager::install("Harman")' && \
    sudo R -e 'BiocManager::install("ImpulseDE2")' && \
    sudo R -e 'BiocManager::install("maSigPro")' && \
    sudo R -e 'BiocManager::install("ENCODExplorer")' && \
    sudo R -e 'BiocManager::install("enrichplot")' && \
    sudo R -e 'BiocManager::install("fgsea")' && \
  #  sudo R -e 'BiocManager::install("IRanges")' && \
  #  sudo R -e 'BiocManager::install("pathview")' && \
  #  sudo R -e 'BiocManager::install("RcisTarget")' && \
  #  sudo R -e 'BiocManager::install("STRINGdb")' && \
  #  sudo R -e 'BiocManager::install("topGO")' && \
    sudo R -e 'install.packages("future")' && \
    sudo R -e 'install.packages("heatmaply")' && \
    sudo R -e 'BiocManager::install("org.Dm.eg.db")' && \
    sudo R -e 'BiocManager::install("org.Hs.eg.db")' && \
    sudo R -e 'BiocManager::install("org.Mm.eg.db")' && \
    sudo R -e 'BiocManager::install("DESeq2")' && \
  #  sudo R -e 'BiocManager::install("BiocParallel")' && \
  #  sudo R -e 'BiocManager::install("biomaRt")' && \   
  #  sudo R -e 'BiocManager::install("DOSE")' && \

    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
    chown shiny:shiny /var/lib/shiny-server
 
# Get pip packages
RUN pip install numpy
RUN pip install pandas
RUN pip install natsort
# RUN pip install multiqc
# RUN pip install pybigwig
# RUN pip install pysam
# RUN pip install pyyaml
# RUN pip install seaborn
# RUN pip install deeptools
# RUN pip install intervene
# RUN pip install htseq # 1) if needed try .tar, then 2) update to python 3 

# Build Bowtie2 (https://sites.google.com/site/wiki4metagenomics/tools/bowtie2/install)
RUN apt update
RUN apt install bowtie2

# Build HISAT2 (https://www.howtoinstall.me/ubuntu/18-04/hisat2/)
RUN apt update
RUN apt install hisat2

# Build samtools
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar -xvzf samtools-1.11.tar.bz2
cd samtools-1.11
./configure --prefix=/where/to/install
make
make install

# Build sra-tools
apt-get --quiet install --yes libxml-libxml-perl
echo "installing sra toolkit to /usr/local/ncbi"
rm -rf .ncbi /usr/local/ncbi /etc/ncbi /etc/profile.d/sra-tools* # remove old install if any
curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-centos_linux64-cloud.tar.gz | tar xz -C /
echo "Please 'source /etc/profile.d/sra-tools.sh' to setup your path"
RUN 'source /etc/profile.d/sra-tools.sh'

# Build fastqc
ENV URL=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV ZIP=fastqc_v0.11.8.zip
RUN wget $URL/$ZIP -O $DST/$ZIP && \
  unzip - $DST/$ZIP -d $DST && \
  rm $DST/$ZIP && \
  cd $DST/FastQC && \
  chmod 755 fastqc && \
  ln -s $DST/FastQC/fastqc /usr/local/bin/fastqc
ENV PATH /usr/local/bin:$PATH

EXPOSE 3838

COPY shiny-server.sh /usr/bin/shiny-server.sh
RUN mkdir -p /srv/demos/
ADD app /srv/shiny-server/
ADD demos /srv/demos/
RUN chown -R shiny:shiny /srv/demos/

CMD ["/usr/bin/shiny-server.sh"]
