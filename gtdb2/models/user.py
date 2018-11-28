from django.db import models
from .abstract import AbstractUnit, AbstractParam


class User(AbstractUnit):
    password = models.CharField(max_length=255, default=None, blank=True, null=True)
    email = models.CharField(max_length=255, default=None, blank=True, null=True)

    class Meta:
        db_table = 'users'


class UserParam(AbstractParam):
    user = models.ForeignKey(User, on_delete=models.CASCADE)

    class Meta:
        db_table = 'user_params'

