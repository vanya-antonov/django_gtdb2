from django.db import models


class User(models.Model):
    c_date = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=255, unique=True)
    descr = models.CharField(max_length=255, default=None, blank=True, null=True)
    password = models.CharField(max_length=255, default=None, blank=True, null=True)
    email = models.CharField(max_length=255, default=None, blank=True, null=True)

    class Meta:
        db_table = 'users'

