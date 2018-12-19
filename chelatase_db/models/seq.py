# Copyright 2018 by Ivan Antonov. All rights reserved.

from gtdb2.models.seq import Seq


class ChelataseSeq(Seq):
    class Meta:
        proxy = True

    # Merge two dicts: https://stackoverflow.com/a/38990/310453
    # prm_info = dict(list(Seq.prm_info.items()) + list({}.items()))

