FROM python:2.7

WORKDIR /code

RUN pip install --upgrade pip
RUN pip install invoke

COPY tasks.py /code

ENTRYPOINT ["invoke"]
