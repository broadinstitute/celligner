FROM python:3.8-alpine
COPY . /app
WORKDIR /app
RUN pip install .
