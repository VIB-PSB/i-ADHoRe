FROM ubuntu

MAINTAINER bedro

COPY . /i-adhore
WORKDIR /i-adhore
RUN apt-get update && apt-get -y install curl && apt-get -y install cmake 
RUN cd /i-adhore
RUN mkdir build 
RUN cd build 
RUN ls


