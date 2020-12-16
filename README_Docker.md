# Follow these steps to build TIMEOR Docker Image
## Ashley Mae Conard
### Last Mod. Aug. 9, 2020

### To build Docker image:

docker build -t timeor_env .

### To run Docker image:

docker run -p 9111:3838 timeor_env

### In this example, R Shiny will be running on port 9111
