#!/bin/sh
# Resource: https://github.com/rocker-org/shiny/blob/master/shiny-server.sh

# Make sure the directory for individual app logs exists
mkdir -p /var/log/shiny-server
chown shiny.shiny /var/log/shiny-server

if [ "$APPLICATION_LOGS_TO_STDOUT" != "false" ];
then
    # push the "real" application logs to stdout with xtail in detached mode
    exec xtail /var/log/shiny-server/ &
fi

head -n 2 /home/shiny/.ncbi/vdb-user-settings.mkfg  > /home/shiny/.ncbi/user-settings.mkfg 
printf '/LIBS/GUID = "%s"\n' `uuidgen` >>  /home/shiny/.ncbi/user-settings.mkfg 
tail -n +4   /home/shiny/.ncbi/vdb-user-settings.mkfg  >>   /home/shiny/.ncbi/user-settings.mkfg 
chmod 600  /home/shiny/.ncbi/user-settings.mkfg  
chown -R shiny:shiny /home/shiny/.ncbi 

# start shiny server
exec shiny-server 2>&1
