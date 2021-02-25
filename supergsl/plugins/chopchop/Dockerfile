FROM python:2.7

WORKDIR /code

COPY requirements.txt /code

RUN pip install --upgrade pip
RUN pip install -r requirements.txt


ADD https://bitbucket.org/valenlab/chopchop/get/bf88f13e387d.zip ./
RUN unzip bf88f13e387d.zip
RUN mv valenlab-chopchop-bf88f13e387d chopchop

COPY tasks.py /code

ENTRYPOINT ["invoke"]
