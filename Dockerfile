FROM rocker/r-ver:3.6.1

RUN apt-get update && apt-get install -y \
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
    wget

# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    sudo R -e 'install.packages(c("openssl", "httr", "rvest", "xml2"), dependencies=TRUE)' && \
    sudo R -e 'install.packages(c("autoplotly","biocmanager","broom","cluster","d3heatmap","data.table","devtools","domc","dorng","dplyr","dt","europepmc","factoextra","farver","fastmatch","fpc","ggforce","ggfortify","ggplot2","ggplotify","ggridges","globaloptions","graphlayouts","gridgraphics","heatmaply","markdown","plotly","png","polyclip","promises","reticulate","rmarkdown","rvcheck","shiny","shinyalert","shinycssloaders","shinydashboard","shinydashboardplus","shinyjs","shinylp","shinyWidgets","stringr","tibble","tidygraph","tidyr","tidyverse","triebeard","tweenr","urltools","vegan","zoo"))' && \
    sudo R -e 'BiocManager::install("Harman")' && \
    sudo R -e 'install.packages("shinydashboardPlus")' && \
    sudo R -e 'install.packages("shinyLP")' && \
    sudo R -e 'install.packages("future")' && \
    sudo R -e 'BiocManager::install("org.Dm.eg.db")' && \
    sudo R -e 'BiocManager::install("DESeq2")' && \ 
    #R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN')" && \
    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
    chown shiny:shiny /var/lib/shiny-server

EXPOSE 3838

#RUN mkdir -p /usr/bin/shiny-server/
COPY shiny-server.sh /Users/ashleyconard/Desktop/shiny-server.sh
ADD app /Users/ashleyconard/Desktop/shiny-server/
ADD demos /Users/ashleyconard/Desktop/shiny-server/

CMD ["/Users/ashleyconard/Desktop/shiny-server.sh"]
