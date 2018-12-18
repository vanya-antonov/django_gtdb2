# Copyright 2018 by Ivan Antonov. All rights reserved.

import json
import logging
import sys

from django.db import models
from django.utils import timezone

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.user import User


class AbstractUnit(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    c_date = models.DateTimeField(default=timezone.now)
    name = models.CharField(max_length=255, unique=True)
    descr = models.CharField(max_length=255, default=None, blank=True, null=True)

    class Meta:
        abstract = True

    gtdb = GeneTackDB()

    prm_info = {
        'xref': {'value_attr': 'data', 'type_fun': json.loads, 'is_list': True},
    }

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
                info = {'value_attr': 'value'}

            # Sort the list by one of the attributes if needed
            if 'sort_attr' in info:
                param_list = sorted(param_list,
                                    key=lambda p: getattr(p, info['sort_attr']),
                                    reverse=info.get('reverse', False))

            # Select one of the param attributes only
            param_list = [getattr(p, info['value_attr']) for p in param_list]

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

    def add_gbk_xref_param(self, gbk_xref):
        "gbk_xref is a string like 'Assembly:GCF_001613165.1'."
        parts = gbk_xref.split(':')
        if len(parts) != 2:
            raise ValueError("Wrong gbk_xref='%s'" % gbk_xref)
        return self.add_xref_param(parts[0], parts[1])

    def add_xref_param(self, db_name, ext_id, data_dict=None):
        "Creates a new param with name='xref', value and data (json)."
        if data_dict is None:
            data_dict = {'db_name': db_name, 'ext_id': ext_id}
        else:
            # Validate the content of data_dict
            if data_dict['db_name'] != db_name:
                raise ValueError("Wrong db_name in data_dict: %s", data_dict)
            if data_dict['ext_id'] != ext_id:
                raise ValueError("Wrong ext_id in data_dict: %s", data_dict)
        json_str = json.dumps(data_dict)
        value = db_name + ':' + str(ext_id)
        return self.add_param('xref', value, data=json_str)


class AbstractParam(models.Model):
    name = models.CharField(max_length=255)
    value = models.CharField(max_length=255, default=None, blank=True, null=True)
    num = models.FloatField(null=True, blank=True, default=None)
    data = models.TextField(null=True, blank=True, default=None)

    class Meta:
        abstract = True

    def __str__(self):
        return self.name

    @classmethod
    def get_parent_by_xref(cls, db_name, ext_id):
        "Returns a Unit object (or None) by external database ID."
        value = db_name + ':' + str(ext_id)
        all_params = cls.objects.filter(name='xref', value=value).all()
        if len(all_params) == 0:
            return None
        elif len(all_params) == 1:
            return all_params[0].parent
        else:
            raise ValueError("xref value '%s' is not unique in model '%s'" %
                             (value, cls.__name__))

