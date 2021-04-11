FROM python:3.8

RUN apt-get update && apt-get install -y \
    docker.io

RUN pip install -U pip

RUN mkdir /code
WORKDIR /code
COPY . /code

RUN python setup.py develop

ENTRYPOINT [ "invoke" ]
