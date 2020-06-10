FROM ubuntu

MAINTAINER bedro

ENV DEBIAN_FRONTEND=noninteractive

COPY . /i-adhore

RUN apt-get update && apt-get install -y \
    curl \
    cmake \
    build-essential

RUN cd /i-adhore \
    && mkdir build 

WORKDIR /i-adhore/build 

RUN cmake .. 
RUN make
RUN make install




