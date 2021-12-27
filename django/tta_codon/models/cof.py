from gtdb2.models import Cof


class TtaCof(Cof):
    class Meta:
        proxy = True

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(Cof.PRM_INFO.items()) + list({
        'tta__num_tta_genes': {'value_attr': 'num', 'type_fun': int},
        'tta__num_tta_fs_genes': {'value_attr': 'num', 'type_fun': int},
    }.items()))

