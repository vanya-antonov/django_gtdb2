# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint
import os

from django.conf import settings

from chelatase_db.models.fshift import ChelataseFshift
#from chelatase_db.models.org import ChelataseOrg
from chelatase_db.models.seq import ChelataseSeq
from gtdb1.models import GtFs
from gtdb2.models.cof import Cof
from gtdb2.models.fshift import FshiftParam
from gtdb2.models.org import Org


class ChelataseCof(Cof):
    class Meta:
        proxy = True

    chld_info = {
        'name': 'chlD',
        'descr': 'Mg/Co-chelatase medium subunit',
        'seed_fshifts': {
            286671156: {'org_gbk': 'chelatase_db/data/D_acidovorans.gbk'},
            33540948: {'org_gbk': 'chelatase_db/data/M_fervens.gbk'},
            946577251: {'org_gbk': 'chelatase_db/data/M_sp.gbk'}}
    }

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    # prm_info = dict(list(Cof.prm_info.items()) + list({
    # }.items()))

    @classmethod
    def get_or_create_chld_cof(cls, user):
        "Returns the chlD cof."
        cof = cls.objects.filter(name=cls.chld_info['name']).first()
        if cof is None:
            gtdb1_fs_list = []
            for gtdb1_id in cls.chld_info['seed_fshifts'].keys():
                gtdb1_fs = GtFs.objects.using('gtdb1').get(pk=gtdb1_id)
                gtdb1_fs_list.append(gtdb1_fs)
            return cls.create_from_gtdb1_seed_fshifts(user, gtdb1_fs_list)
        else:
            if cof.descr != cls.chld_info['descr']:
                raise ValueError("Name/descr mismatch!")
            return cof

    @classmethod
    def create_from_gtdb1_seed_fshifts(cls, user, gtdb1_fshifts):
        """Creates the chlD cof with the seed frameshifts from GTDB1.

        Arguments:
         - gtdb1_fshifts - list of GtFs objects loaded from GTDB1.
        """
        seed_fshifts = []
        for gtdb1_fs in gtdb1_fshifts:
            fshift = FshiftParam.get_parent_by_xref('gtdb1', gtdb1_fs.fs_id)
            if fshift is None:
                info = cls.chld_info['seed_fshifts'][gtdb1_fs.fs_id]
                gbk_fn = os.path.join(settings.BASE_DIR, info['org_gbk'])

                # Use Org instead of ChelataseOrg so that it doesn't try to
                # find chld genes using chld_cof that doesn't exist yet
                #org = ChelataseOrg.get_or_create_from_gbk(user, gbk_fn)
                org = Org.get_or_create_from_gbk(user, gbk_fn)

                seq = ChelataseSeq.get_or_create_from_ext_id(
                    user, org, gtdb1_fs.job.name + '.1')
                seq.make_all_params()

                fshift = ChelataseFshift.create_from_gtdb1_fs(
                    user, seq, gtdb1_fs)
                fshift.make_all_params()

            #if fshift.feat is None:
            #    ChelataseFeat.create_fsCDS_from_fshift(user, fshift)

            seed_fshifts.append(fshift)

        return cls.create_from_seed_fshifts(
            user, seed_fshifts, name=cls.chld_info['name'],
            descr=cls.chld_info['descr'])

