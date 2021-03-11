FROM python:2.7

WORKDIR /code

COPY requirements.txt /code

RUN pip install --upgrade pip

COPY tasks.py /code

ENTRYPOINT ["invoke"]
