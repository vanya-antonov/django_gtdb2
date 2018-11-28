from django.db import models


class AbstractUnit(models.Model):
    c_date = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=255, unique=True)
    descr = models.CharField(max_length=255, default=None, blank=True, null=True)

    class Meta:
        abstract = True

    def __str__(self):
        return self.name


class AbstractParam(models.Model):
    name = models.CharField(max_length=255)
    value = models.CharField(max_length=255, default=None, blank=True, null=True)
    num = models.FloatField(null=True, blank=True, default=None)
    data = models.TextField(null=True, blank=True, default=None)

    class Meta:
        abstract = True

