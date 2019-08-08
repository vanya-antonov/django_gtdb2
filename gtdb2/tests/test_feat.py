# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint

from gtdb2.models.org import Org
from gtdb2.models.seq import Seq
from gtdb2.models.feat import Feat
from gtdb2.models.fshift import Fshift
from gtdb2.tests import GtdbTestCase


class FeatModelTests(GtdbTestCase):

    def setUp(self):
        super().setUp()
        fn = self.get_full_path_to_test_file('M_fervens.gbk')
        self.org = Org.create_from_gbk(self.user, fn)
        self.seq = Seq.create_from_ext_id(
            self.user, self.org, 'NC_013156.1')

    def test_feat_create_from_gbk_annotation(self):
        start, end = 1121860-1, 1122960   # convert to 0-based coordinates
        feat = Feat.get_or_create_from_gbk_annotation(
            self.user, self.seq, start, end, 1)

        # Check attributes
        self.assertEqual(feat.name, 'MEFER_RS06095')
        self.assertEqual(feat.descr, 'magnesium chelatase')
        self.assertEqual(feat.type, 'CDS')
        self.assertEqual(feat.origin, 'annotation')

        # Check params
        self.assertEqual(feat.prm.location_str, 'NC_013156.1:1121859-1122960(+)')
        self.assertEqual(feat.prm['protein_id'], 'WP_015791742.1')
        self.assertTrue(feat.prm['translation'].startswith('MQYIYPFTAIVGQ'))
        # https://www.ncbi.nlm.nih.gov/nuccore/NC_013156.1?report=fasta&from=1121860&to=1122960
        self.assertTrue(feat.prm.seq_nt.startswith('ATGCAATAT'))
        self.assertTrue(feat.prm.seq_nt.endswith('AATTAA'))

        # Make sure new feature is not created if we call the function again
        # and the existing Feat is returned
        self.assertEqual(Feat.objects.filter(type='CDS').count(), 1)
        feat2 = Feat.get_or_create_from_gbk_annotation(
            self.user, self.seq, start+10, end-20, 1)
        self.assertEqual(Feat.objects.filter(type='CDS').count(), 1)
        self.assertEqual(feat, feat2)

    def test_feat_get_or_create_fscds_from_parent(self):
        start, end = 1121860-1, 1122960   # convert to 0-based coordinates
        parent_feat = Feat.get_or_create_from_gbk_annotation(
            self.user, self.seq, start, end, 1)

        # fsCDS needs a frameshift for full length translation
        fshift = Fshift.get_or_create(
            user=self.user, seq=self.seq, origin='tblastn',
            start=1121859, end=1123904, strand=1,
            coord=1122957, len=-1)

        # Finally, create the full-length fsCDS feat
        fscds = Feat.get_or_create_fscds_from_parent(
            self.user, parent_feat, [fshift])

        self.assertEqual(fscds.type, 'fsCDS')
        self.assertEqual(fscds.origin, 'genetack')
        self.assertEqual(fscds.name, 'MEFER_RS06095_fs1122957')
        self.assertEqual(list(fscds.fshift_set.all()), [fshift])

        # Check the fsCDS params
        self.assertTrue(fscds.prm.translation.startswith('MQYIYPFTA'))
        self.assertTrue(fscds.prm.translation.endswith('ICKGFVEY'))

        # https://www.ncbi.nlm.nih.gov/nuccore/NC_013156.1?report=fasta&from=1121859&to=1123904
        self.assertTrue(fscds.prm.seq_nt.startswith('ATGCAATAT'))
        self.assertTrue(fscds.prm.seq_nt.endswith('TATTAG'))

        # Make sure new feature is not created if we call the function again
        # and the existing Feat is returned
        self.assertEqual(Feat.objects.filter(type='fsCDS').count(), 1)
        fscds2 = Feat.get_or_create_fscds_from_parent(
            self.user, parent_feat, [fshift])
        self.assertEqual(Feat.objects.filter(type='fsCDS').count(), 1)
        self.assertEqual(fscds, fscds2)

