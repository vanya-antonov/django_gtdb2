FROM python:3.7
RUN mkdir /app_dir
WORKDIR /app_dir

ENV PYTHONUNBUFFERED 1

COPY ./requirements.txt /app_dir
RUN pip install --upgrade pip && pip install -r requirements.txt
COPY . /app_dir
