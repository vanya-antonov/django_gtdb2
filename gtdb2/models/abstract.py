import logging
import sys

from django.db import models

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.user import User


class AbstractUnit(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    c_date = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=255, unique=True)
    descr = models.CharField(max_length=255, default=None, blank=True, null=True)

    class Meta:
        abstract = True

    gtdb = GeneTackDB()
    prm_info = {}

    def __str__(self):
        return self.name

    def _get_param_set(self):
        # e.g. unit = 'org'
        unit = self.__class__.__name__.lower()
        return getattr(self, unit + "param_set")
    param_set = property(
        fget=_get_param_set,
        doc="Convenience property: e.g. org.param_set == org.orgparam_set")

    def _get_param_dict(self):
        param_dict = {}
        for p in self.param_set.all():
            param_dict.setdefault(p.name, []).append(p)
        return param_dict
    param_dict = property(
        fget=_get_param_dict,
        doc="""Convenience property:
        >>> list_A = list(org.orgparam_set.filter(name='source_fn').all())
        >>> list_B = org.param_dict['source_fn']
        >>> list_A == list_B
        """)

    def _get_prm(self):
        prm = {}
        for key, param_list in self.param_dict.items():
            info = self.prm_info.get(key, None)
            if info is None:
                continue

            # Sort the list by one of the attributes if needed
            if 'sort_attr' in info:
                param_list = sorted(param_list,
                                    key=lambda p: getattr(p, info['sort_attr']),
                                    reverse=info.get('reverse', False))

            # Select one of the param attributes only
            param_list = [getattr(p, info['prm_attr']) for p in param_list]

            # Convert it to another type if needed
            if 'type_fun' in info:
                type_fun = info['type_fun']
                param_list = [type_fun(p) for p in param_list]

            if info.get('is_list', False):
                prm[key] = param_list
            else:
                if len(param_list) > 1:
                    logging.error("param_list contains '%s' elements for the "
                                  "non list prm '%s'" % (len(param_list)), key)
                prm[key] = param_list[0]
        return prm
    prm = property(
        fget=_get_prm,
        doc="A computed dictionary with parameters represented by simple types.")

    def set_param(self, name, value=None, num=None, data=None):
        "Deletes and adds param."
        self.delete_param(name)
        self.add_param(name, value, num, data)

    def add_param(self, name, value=None, num=None, data=None):
        """E.g. for Org object it will call:
        >>> param = OrgParam(...)
        >>> param.save()
        """
        # https://stackoverflow.com/a/41236263/310453
        # Org => 'gtdb2.models.org'
        modulename = "gtdb2.models." + self.__class__.__name__.lower()
        mod = sys.modules[modulename]

        # Org => 'OrgParam'
        classname = self.__class__.__name__ + 'Param'
        cls = getattr(mod, classname)

        param = cls(parent=self, name=name, value=value, num=num, data=data)
        param.save()
        return param

    def delete_param(self, name):
        "Deletes all unit params with the given name."
        self.param_set.filter(name=name).delete()


class AbstractParam(models.Model):
    name = models.CharField(max_length=255)
    value = models.CharField(max_length=255, default=None, blank=True, null=True)
    num = models.FloatField(null=True, blank=True, default=None)
    data = models.TextField(null=True, blank=True, default=None)

    class Meta:
        abstract = True

    def __str__(self):
        return self.name

