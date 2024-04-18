FROM ubuntu:20.04

MAINTAINER bedro

ENV DEBIAN_FRONTEND=noninteractive

COPY . /i-adhore

RUN apt-get update && apt-get install -y \
    curl \
    cmake \
    libpng-dev \
    build-essential

RUN cd /i-adhore \
    && mkdir build 

WORKDIR /i-adhore/build 

RUN cmake .. 
RUN make
RUN make install

WORKDIR /

