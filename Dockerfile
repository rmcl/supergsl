FROM python:3.9

RUN apt-get update && apt-get install -y \
    docker.io

RUN pip install -U pip

RUN mkdir /code
WORKDIR /code
COPY . /code

RUN pip install -e .[plugins]

RUN python setup.py develop

ENTRYPOINT [ "invoke" ]
