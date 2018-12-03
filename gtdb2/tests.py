import os

from django.conf import settings
from django.test import TestCase

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.org import Org


class OrgModelTests(TestCase):
    def test_org_create_from_gbk(self):
        db = GeneTackDB()
        user = db.get_default_user()
        fn = os.path.join(settings.BASE_DIR, 'gtdb2/data/N_mexicana.gbk')

        org = Org.create_from_gbk(user, fn)

        # org.name is species name, so it must include genus, e.g.
        # org.genus = 'Nocardia'
        # org.name = 'Nocardia mexicana NBRC 108244'
        self.assertIn(org.genus, org.name)

        # Make sure the fn was saved as one of the params
        self.assertEqual(org.orgparam_set.get(name='source_fn').value, fn)

        # Test AbstractUnit computed properties
        self.assertEqual(org.param_set.get(name='source_fn').value, fn)
        self.assertEqual(org.param_dict['source_fn'][0].value, fn)
        #self.assertEqual(org.prm['source_fn'], fn)


