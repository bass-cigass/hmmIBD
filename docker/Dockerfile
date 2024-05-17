FROM alpine as build-env
RUN apk add --no-cache build-base
RUN apk add --no-cache --upgrade bash
RUN apk add --no-cache python3 py3-pip

RUN apk add --no-cache python3-dev

RUN pip install matplotlib --break-system-packages
RUN pip install numpy --break-system-packages
RUN pip install seaborn --break-system-packages
copy py py 

Copy ./hmmIBD hmmIBD
WORKDIR ./hmmIBD

# Compile the binaries
RUN gcc -o hmmIBD hmmIBD.c

RUN cp hmmIBD /usr/bin/