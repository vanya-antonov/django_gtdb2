FROM node:latest as build-stage
ARG VUE_APP_API_BASE_URL=$VUE_APP_API_BASE_URL
WORKDIR /app
COPY frontend/package*.json ./
RUN npm install
COPY ./frontend/ .
RUN npm run build  --loglevel verbose

FROM nginx:1.17.10 as production-stage
RUN mkdir /app
RUN mkdir /app/chelatase
COPY --from=build-stage /app/dist /app
