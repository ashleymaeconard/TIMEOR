FROM jlaw2677/timeor_most_recent
RUN wget https://github.com/DaehwanKimLab/hisat2/archive/0f01dc6397a.tar.gz && tar -xvzf 0f01dc6397a.tar.gz && cd hisat2-0f01dc6397a && make
RUN mkdir -p /srv/demos/
ADD app /srv/shiny-server/
ADD demos /srv/demos/
RUN chown -R shiny:shiny /srv/demos/
RUN apt-get update && apt-get install uuid-runtime
RUN su - shiny && mkdir -p /home/shiny/.ncbi
COPY vdb-user-settings.mkfg /home/shiny/.ncbi
COPY shiny-server_2.sh /usr/bin/shiny-server_2.sh 
# CMD  head -n 2 /home/shiny/.ncbi/vdb-user-settings.mkfg  > /home/shiny/.ncbi/user-settings.mkfg && \
# printf '/LIBS/GUID = "%s"\n' `uuidgen` >>  /home/shiny/.ncbi/user-settings.mkfg && \
# tail -n +4   /home/shiny/.ncbi/vdb-user-settings.mkfg  >>   /home/shiny/.ncbi/user-settings.mkfg && \
# chmod 600  /home/shiny/.ncbi/user-settings.mkfg  && chown shiny:shiny  /home/shiny/.ncbi/user-settings.mkfg && \
CMD ["/usr/bin/shiny-server_2.sh"]
