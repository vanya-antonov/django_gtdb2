FROM node:latest as build-stage
WORKDIR /app/chelatase
COPY frontend/package*.json ./
RUN npm install
COPY ./frontend/ .
RUN npm run build

FROM nginx:1.17.10 as production-stage
RUN mkdir /app
RUN mkdir /app/chelatase
COPY --from=build-stage /app/chelatase/dist /app/chelatase
