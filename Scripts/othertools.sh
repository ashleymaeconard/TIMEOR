
# Build Bowtie2 (https://sites.google.com/site/wiki4metagenomics/tools/bowtie2/install)
apt update
apt install bowtie2

# Build HISAT2 (https://www.howtoinstall.me/ubuntu/18-04/hisat2/)
apt update
apt install hisat2

# Build samtools 
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar xvjf samtools-1.11.tar.bz2 && \
    rm  samtools-1.11.tar.bz2 && \
    cd samtools-1.11 && \
    ./configure --without-curses && \
    make && \
    make install && \
    cd ..

# Build sra-tools
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -vxzf sratoolkit.tar.gz && \
    PATH=$PATH:$PWD/sratoolkit.2.10.9-ubuntu64/bin && \
    rm sratoolkit.tar.gz 

# Build fastqc
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip \
    && unzip fastqc_v0.11.8.zip \
    && rm fastqc_v0.11.8.zip \
    && chmod 755 /FastQC/fastqc 

# Build deeptools
wget https://github.com/deeptools/deepTools/archive/1.5.12.tar.gz && \
    tar -xzvf 1.5.12.tar.gz && \
    rm 1.5.12.tar.gz && \
    cd /deepTools-1.5.12 && \
    python setup.py install && \
    cd ..

# Build htseq
wget --no-check-certificate https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1p1.tar.gz && \
    tar -zxvf HTSeq-0.6.1p1.tar.gz && \
    rm HTSeq-0.6.1p1.tar.gz && \
    cd HTSeq-0.6.1p1/ && \
    python setup.py install && \
    chmod +x scripts/htseq-count && \
    chmod +x scripts/htseq-qa
