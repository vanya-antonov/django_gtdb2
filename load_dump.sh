if [ -z "$SQL_DUMP_PATH" ]
then
	echo SET SQL_DUMP_PATH variable!
else
	sudo docker-compose exec -T db sh -c 'mysql -u$MYSQL_USER -p$MYSQL_PASSWORD $MYSQL_DATABASE' <  $SQL_DUMP_PATH
fi
