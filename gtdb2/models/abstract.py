# Copyright 2018 by Ivan Antonov. All rights reserved.

import json
import logging
import sys

from attrdict import AttrDict
from django.db import models
from django.utils import timezone

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.user import User


class AbstractUnit(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    c_date = models.DateTimeField(default=timezone.now)
    name = models.CharField(max_length=255, default=None, blank=True, null=True)
    descr = models.CharField(max_length=255, default=None, blank=True, null=True)

    class Meta:
        abstract = True

    gtdb = GeneTackDB()

    # All items in the 'prm' dict must be explicitly declared in 'prm_info'
    prm_info = {
        'xref': {'value_attr': 'data', 'type_fun': json.loads, 'is_list': True},
    }

    def __str__(self):
        return str(self.name)

    @property
    def param_dict(self):
        """A convenience property:
        >>> list_A = list(org.orgparam_set.filter(name='source_fn').all())
        >>> list_B = org.param_dict['source_fn']
        >>> list_A == list_B
        """
        param_dict = {}
        for p in self.param_set.all():
            param_dict.setdefault(p.name, []).append(p)
        return param_dict

    @property
    def prm(self):
        "A computed dictionary with parameters represented by simple types."
        prm = AttrDict()
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
            value_attr = info.get('value_attr', 'value')
            value_list = [getattr(p, value_attr) for p in param_list]

            # Convert it to another type if needed
            if 'type_fun' in info:
                type_fun = info['type_fun']
                value_list = [type_fun(p) for p in value_list]

            if info.get('is_list', False):
                prm[key] = value_list
            else:
                if len(value_list) > 1:
                    logging.error("value_list contains '%s' elements for the "
                                  "non-list prm '%s'" % (len(value_list)), key)
                prm[key] = value_list[0]
        return prm

    def set_param(self, name, value=None, num=None, data=None):
        "Deletes and adds param."
        self.delete_param(name)
        self.add_param(name, value, num, data)

    def add_param(self, name, value=None, num=None, data=None):
        """E.g. for Org object it will call:
        >>> param = OrgParam(...)
        >>> param.save()
        """
        # E.g. self.param_set.model is 'OrgParam' class:
        # https://stackoverflow.com/a/14487781/310453
        param_cls = self.param_set.model

        # Do not add extra rows to non-list params
        info = self.prm_info.get(name, None)
        if(info is not None and not info.get('is_list', False) and
          param_cls.objects.filter(parent=self, name=name).count() > 0):
            raise ValueError("Can't add another row to a non-list param '%s' "
                             "(parent_id = '%s')" % (name, self.id))

        param = self.param_set.model(
            parent=self, name=name, value=value, num=num, data=data)
        param.save()

        return param

    def delete_param(self, name):
        "Deletes all unit params with the given name."
        self.param_set.filter(name=name).delete()

    def set_param_xref(self, db_name, ext_id):
        "Deletes xref(s) for the given db_name and adds the new ext_id."
        self.param_set.filter(
            name='xref', value__startswith=db_name + ':'
        ).delete()
        self.add_param_xref(db_name, ext_id)

    def add_param_xref(self, db_name, ext_id):
        "Creates a new param with name='xref' and value='db_name:ext_id'."
        data_dict = {'db_name': db_name, 'ext_id': ext_id}
        json_str = json.dumps(data_dict)
        value = db_name + ':' + str(ext_id)
        return self.add_param('xref', value, data=json_str)

    # TODO
    # def add_param_xref_data_dict(self, data_dict):


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

