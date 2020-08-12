FROM continuumio/miniconda
COPY packages.yml /root/
RUN conda env create -f /root/packages.yml
RUN conda init bash
# Activate the environment, and make sure it's activated:
RUN echo "conda activate timeor_conda_env" > ~/.bashrc
RUN conda run -n timeor_conda_env R -e "BiocManager::install('clusterProfiler')"
RUN mkdir -p /var/lib/shiny-server/bookmarks/shiny
# Copy the app to the image
COPY app /root/app
#COPY TIMEOR_tutorial /root/
COPY demos /root/
COPY run.sh /root/
#COPY Rprofile.site /usr/local/lib/R/etc/Rprofile.site
EXPOSE 3838
#CMD ["conda", "run", "-n", "timeor_env_basics", "/bin/bash", "/root/run.sh"]
#CMD ["conda", "run", "-n", "timeor_env_basics", "R", "-e", "print(c('hello', getwd()))"]
CMD ["conda", "run", "-n", "timeor_conda_env", "R", "-e", "shiny::runApp('/root/app', launch.browser=F, host='0.0.0.0', port=3838)"]
#CMD ["/bin/bash"]#conda", "run", "-n", "timeor_env_basics", "R"]
#ENTRYPOINT ["conda", "run", "-n", "timeor_env_basics", "conda", "info"]

