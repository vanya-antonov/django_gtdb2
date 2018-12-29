# Copyright 2018 by Ivan Antonov. All rights reserved.

from gtdb2.models.fshift import Fshift


class ChelataseFshift(Fshift):
    class Meta:
        proxy = True

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    # prm_info = dict(list(Fshift.prm_info.items()) + list({
    # }.items()))

    @property
    def feat(self):
        "Simulate one-to-one relationship between chelatase objects."
        num = self.feat_set.count()
        if num == 0:
            return None
        elif num == 1:
            return self.feat_set.first()
        else:
            raise ValueError("ChelataseFshift '%s' has '%s' feats" %
                             (self, num))

