# Follow these steps to build and/or run TIMEOR Docker image
## Ashley Mae Conard
#### Last Mod. Aug. 9, 2020

### Assuming you have Docker Hub, to pull down Docker image:
`$docker pull ashleymaeconard/timeor:latest`

### Assuming you have Docker installed (version 20.10.0 recommended), to build Docker image:

`$docker build -t timeor_env .`

### To run the Docker image:

`$docker run -p 9111:3838 timeor_env`

### In this example, R Shiny will be running on port 9111
