# Copyright 2018 by Ivan Antonov. All rights reserved.

from chelatase_db.models.fshift import ChelataseFshift
from gtdb1.models import GtFs
from gtdb2.models.cof import Cof
from gtdb2.models.fshift import FshiftParam


class ChelataseCof(Cof):
    class Meta:
        proxy = True

    info = {
        'name': 'chlD',
        'descr': 'Mg/Co-chelatase medium subunit',
        'seed_fshifts': [
            {'gtdb1_id': 286671156, 'org_gbk': 'chelatase_db/data/D_acidovorans.gbk'},
            {'gtdb1_id': 33540948, 'org_gbk': 'chelatase_db/data/M_fervens.gbk'},
            {'gtdb1_id': 946577251, 'org_gbk': 'chelatase_db/data/M_sp.gbk'}]
    }

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    # prm_info = dict(list(.prm_info.items()) + list({
    # }.items()))

    # Overwrite the computed property to keep using the parent's param table
    param_set = property(fget=lambda self: self.cofparam_set)

    @classmethod
    def get_or_create(cls, user):
        cof = cls.objects.filter(name=cls.info['name']).first()
        if cof is not None:
            if cof.descr != cls.info['descr']:
                raise ValueError("Name/descr mismatch!")
            return cof

        seed_fshifts = []
        for fs_info in cls.info['seed_fshifts']:
            fshift = FshiftParam.get_parent_by_xref(
                'gtdb1', fs_info['gtdb1_id'])
            if fshift is None:
                gbk_fn = os.path.join(settings.GTDB_DIR, info['org_gbk'])
                org = ChelataseOrg.get_or_create_from_gbk(user, gbk_fn)
                seq = ChelataseSeq.get_or_create_from_ext_id(user, org, ext_id)

                fs = GtFs.objects.get(pk=info['gtdb1_id'])
                fshift = ChelataseFshift.create_from_gtdb1_fs(user, seq, fs)
                fshift.make_all_params()

            if fshift.feat is None:
                ChelataseFeat.create_fsCDS_from_fshift(user, fshift)

            seed_fshifts.append(fshift)
        return cls.create_from_seed_fshifts(
            seed_fshifts, name=cls.name, descr=cls.descr)

