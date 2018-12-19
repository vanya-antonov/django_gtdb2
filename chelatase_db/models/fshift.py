# Copyright 2018 by Ivan Antonov. All rights reserved.

from gtdb2.models.fshift import Fshift


class ChelataseFshift(Fshift):
    class Meta:
        proxy = True

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    # prm_info = dict(list(Fshift.prm_info.items()) + list({
    # }.items()))

