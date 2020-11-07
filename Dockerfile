FROM python:3.8

RUN pip install -U pip

RUN mkdir /code
WORKDIR /code
COPY . /code

RUN python setup.py develop

ENTRYPOINT [ "invoke" ]
