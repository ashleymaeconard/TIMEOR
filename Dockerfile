FROM continuumio/miniconda2
COPY packages.yml /root/
RUN conda env create -f /root/packages.yml
RUN conda init bash
RUN echo "conda activate timeor_conda_env" > ~/.bashrc
RUN conda run -n timeor_conda_env R -e "BiocManager::install('clusterProfiler')"
RUN mkdir -p /var/lib/shiny-server/bookmarks/shiny
COPY app /root/app
#COPY demos /root/
EXPOSE 3838
CMD ["conda", "run", "-n", "timeor_conda_env", "R", "-e", "shiny::runApp('/root/app', launch.browser=F, host='0.0.0.0', port=3838)"]
