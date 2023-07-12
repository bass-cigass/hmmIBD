FROM alpine as build-env
RUN apk add --no-cache build-base
RUN apk add --no-cache --upgrade bash
RUN apk add --no-cache python3 py3-pip

copy py py 

Copy ./hmmIBD hmmIBD
WORKDIR ./hmmIBD

# Compile the binaries
RUN gcc -o hmmIBD hmmIBD.c

RUN cp hmmIBD /usr/bin/