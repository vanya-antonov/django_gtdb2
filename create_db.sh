#!/bin/bash
if [ -f production.env ]
then
    echo production.env exists
else
	cat prod.env > production.env
	echo MYSQL_PASSWORD=$(base64 /dev/urandom | head -c20) >> production.env
	echo SECRET_KEY=$(base64 /dev/urandom | head -c50) >> production.env
fi
source ./production.env
sudo docker-compose up -d db
