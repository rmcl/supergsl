FROM python:3.8

RUN pip install -U pip

COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

RUN mkdir /code
WORKDIR /code
COPY . /code

ENTRYPOINT [ "invoke" ]
