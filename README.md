# django_gtdb2
GeneTackDB-2 based on Django and Biopython

# Run Django test server instead of Apache on port 80
```bash
sudo systemctl status apache2
sudo systemctl stop apache2.service
sudo python3  ./manage.py  runserver  0.0.0.0:80
```

# Installation on a new machine
```bash

# Install blast
sudo apt-get install ncbi-blast+

git clone https://vanya-antonov@github.com/vanya-antonov/django_gtdb2
cd django_gtdb2

# scp ivan@10.0.1.254:~/_my/github/django_gtdb2/mysite/local_settings.py  mysite/
# - OR -
# scp -P 22194 ivan@83.149.211.146:~/_my/github/django_gtdb2/mysite/local_settings.py  mysite/

python3 -m venv  venv
source venv/bin/activate

# https://stackoverflow.com/a/44862371/310453
pip install --upgrade pip
pip install wheel

# On Ubuntu for mysqlclient:  https://stackoverflow.com/a/7475296/310453
sudo apt-get install libmysqlclient-dev

# On Ubuntu for numpy: https://stackoverflow.com/a/24892867/310453
sudo apt-get install python3.8-dev

cd django
pip install -r requirements.txt

# Copy and edit the local_settings.py
cp -v /home/gtdb/data/local_settings.py  mysite/

# Run server with external access:
# ./manage.py runserver 0.0.0.0:8000
# Run server at: http://127.0.0.1:8000/chelatase_db/
./manage.py runserver

# Test
./manage.py export_seq  translation  Feat  238533
#manage_gtdb2 export_seq  translation  Feat  29451

deactivate
```

# Transfer MySQL database to a new machine
```bash
# [ssh_yulia] Export MySQL data:
mysqldump -u ivan -p -h 127.0.0.1 chelatase_db | gzip > 200417.chelatase_db.sql.gz

# Copy the data file from ssh_yulia
scp -P 22194 ivan@83.149.211.146:~/_my/backup/200417.chelatase_db.sql.gz  .

# On ubuntu the root login is a bit tricky: https://stackoverflow.com/a/42742610/310453
sudo mysql -u root  -p
create database chelatase_db;
create user 'antonov'@'%' identified by '***';
grant all on chelatase_db.* to 'antonov'@'%';

# test connect
alias chelatase_db='mysql -u antonov -p*** chelatase_db'

# About TIMESTAMP NULL:  https://stackoverflow.com/a/37420863/310453
zcat 200417.chelatase_db.sql.gz  |\
perl -pe 's/DATETIME\S+/TIMESTAMP NULL/gi' |\
chelatase_db
```

# Deploy with docker-compose
```bash
cd django_gtdb2
#create file prod.env with database name and mysql_user:
touch prod.env
echo "MYSQL_DATABASE=chelatase_db" >> prod.env
echo "MYSQL_USER=gtdb" >> prod.env
# production.env with secrets will be created automatically
# create production.env and up db
sudo bash create_db.sh
# wait for mysql is ready
sudo docker-compose exec db sh -c 'mysql -u$MYSQL_USER -p$MYSQL_PASSWORD $MYSQL_DATABASE' # should show mysql shell
# unzip sql.gz file
# gunzip ../dump.sql.gz
# load mysql dump
SQL_DUMP_PATH=../170620.chelatase_db.sql bash load_dump.sh
# check loaded tables
sudo docker-compose exec db sh -c 'mysql -u$MYSQL_USER -p$MYSQL_PASSWORD $MYSQL_DATABASE -e "show tables"'

# Start all the containers & go to http://localhost/chelatase
docker-compose up -d
# Check that they are running:
docker ps -a
# Stop all the containers:
docker-compose down
```

# Run unit tests
```bash
# -k option keeps the test database for faster testing
python3 manage.py test -k gtdb2
python3 manage.py test -k chelatase_db

# One test only
./manage.py test -k gtdb2.tests.test_org.OrgModelTests.test_org_run_genetack_gm
```
