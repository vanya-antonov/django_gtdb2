# django_gtdb2
GeneTackDB-2 based on Django and Biopython

# Installation on a new machine
```bash
git clone https://vanya-antonov@github.com/vanya-antonov/django_gtdb2
cd django_gtdb2

# Copy and edit the local_settings.py
# scp ivan@10.0.1.254:~/_my/github/django_gtdb2/mysite/local_settings.py  mysite/
# - OR -
# scp -P 22194 ivan@83.149.211.146:~/_my/github/django_gtdb2/mysite/local_settings.py  mysite/

python3 -m venv  venv
source venv/bin/activate

# https://stackoverflow.com/a/44862371/310453
pip install wheel
pip install -r requirements.txt

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
