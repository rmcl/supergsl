FROM fsharp:10.10.0

RUN mkdir /app
WORKDIR /app

RUN  apt-get update -y --fix-missing && \
     apt-get upgrade -y && \
     apt-get dist-upgrade -y && \
     apt-get -y autoremove && \
     apt-get clean
RUN apt-get install -y \
    zip \
    patch \
    python3 \
    python3-pip

RUN pip3 install invoke

### INSTALL GSL

RUN echo "Fetching GSLc from Github"
ADD https://github.com/Amyris/Gslc/archive/master.zip master.zip
RUN unzip master.zip && mv "Gslc-master" "Gslc"

RUN echo "Building GSLc"
WORKDIR /app/Gslc

COPY build_fsx.patch .
RUN patch -u -b build.fsx -i build_fsx.patch

RUN VisualStudioVersion=14.0 ./build.sh

WORKDIR /app

COPY tasks.py .

ENTRYPOINT ["invoke"]
