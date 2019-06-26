# django_gtdb2
GeneTackDB-2 based on Django and Biopython

# Installation on a new machine
```bash
git clone https://vanya-antonov@github.com/vanya-antonov/django_gtdb2

# Copy and edit the local_settings.py
scp ivan@10.0.1.254:~/_my/github/django_gtdb2/mysite/local_settings.py  mysite/

python3 -m venv  venv
source venv/bin/activate

# https://stackoverflow.com/a/44862371/310453
pip install wheel
pip install -r requirements.txt

# Test
./manage.py export_seq  translation  feat  29451
manage_gtdb2 export_seq  translation  feat  29451

deactivate
```
