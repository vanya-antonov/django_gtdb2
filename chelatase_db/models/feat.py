# Copyright 2018 by Ivan Antonov. All rights reserved.

from gtdb2.models.feat import Feat


class ChelataseFeat(Feat):
    class Meta:
        proxy = True

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    # prm_info = dict(list(Feat.prm_info.items()) + list({
    # }.items()))

    @property
    def fshift(self):
        "Simulate one-to-one relationship between chelatase objects."
        num = self.fshift_set.count()
        if num == 0:
            return None
        elif num == 1:
            return self.fshift_set.first()
        else:
            raise ValueError("ChelataseFeat '%s' has '%s' fshifts" %
                             (self, num))

    @classmethod
    def get_or_create_from_fshift(cls, user, fshift):
        "The fshift should be ChelataseFshift object."
        raise NotImplementedError("Feat!")
        if fshift.feat is not None:
            return fshift.feat
        else:
            return cls.create_from_fshift(user, fshift)

