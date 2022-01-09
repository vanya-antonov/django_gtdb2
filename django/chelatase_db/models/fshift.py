# Copyright 2018 by Ivan Antonov. All rights reserved.
import re

import numpy as np

from pprint import pprint

from chelatase_db.lib.bio import get_stop_stop_seq
from gtdb2.models.fshift import Fshift


def _check_stop_codon_after_fshift(
    seq, strand, parent_feat_coord, fshift_coord, fshift_len, transl_table
):
    feat_start, feat_end = parent_feat_coord
    if strand == 1:
        nearest_stop_ = get_stop_stop_seq(
            feat_end + fshift_len - 3, feat_end + fshift_len, strand, seq, transl_table
        )
        if nearest_stop_["left"] is not None and nearest_stop_["left"] > fshift_coord:
            return True
    if strand == -1:
        nearest_stop_ = get_stop_stop_seq(
            feat_start - fshift_len, feat_start - fshift_len + 3, strand, seq, transl_table,
        )
        if nearest_stop_["right"] is not None and nearest_stop_["right"] < fshift_coord:
            return True
    return False


def search_poly_a_slippery(seq, strand, parent_feat_coord, fshift_len, transl_table):
    start, end = parent_feat_coord
    seq = seq.lower()
    if strand == -1:
        pattern = re.compile(r"t{5,20}")
        match = pattern.search(str(seq), start, end)
        if not match:
            return None
        span = match.span()
        if span[1] > end:
            return None
        fshift_coord = span[1]
        assert set(seq[slice(*span)]) == {"t"}
        assert seq[span[1]] != "t"
        assert seq[span[0] - 1] != "t"

    if strand == 1:
        pattern = re.compile(r"a{5,20}")
        match = pattern.search(str(seq)[::-1], len(seq) - end)
        if not match:
            return None
        span = match.span()
        span = (len(seq) - span[1], len(seq) - span[0])
        if span[0] < start:
            return None
        assert set(seq[slice(*span)]) == {"a"}
        assert seq[span[1]] != "a"
        assert seq[span[0] - 1] != "a"
        fshift_coord = span[0]
        check = _check_stop_codon_after_fshift(
            seq, strand, parent_feat_coord, fshift_coord, fshift_len, transl_table
        )
        if check:
            return None
    return span


def search_stable_secondary_structure(seq, min_length, max_length):
    from rna_tools.Seq import RNASequence

    lengths = []
    energies = []
    preds = []
    for lenth in range(30, 201, 1):
        rna_seq = RNASequence(seq[:lenth])
        pred = rna_seq.predict_ss()
        line, energy = pred.split(" ", 1)
        lengths.append(lenth)
        preds.append(line)
        energies.append(float(energy[1:-1]))
    lenths = np.array(lengths)
    energies = np.array(energies)
    preds = np.array(preds)
    norm_energies = energies / lenths
    ind = norm_energies.argmin()
    pred = preds[ind]
    energy = norm_energies[ind]
    lenth = lenths[ind]
    return {
        "pred": pred,
        "length": lenth,
        "seq": seq[:lenth],
        "energy": energy,
        "preds": preds,
        "lengths": lenths,
        "energies": energies,
    }


class ChelataseFshift(Fshift):
    class Meta:
        proxy = True

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    PRM_INFO = {
        **Fshift.PRM_INFO,
        "poly_a_start": {"value_attr": "num", "type_fun": int},
        "poly_a_end": {"value_attr": "num", "type_fun": int},
        "signal_ss_start": {"value_attr": "num", "type_fun": int},
        "signal_ss_end": {"value_attr": "num", "type_fun": int},
        "signal_ss_struct": {"value_attr": "data"},
        "signal_ss_seq": {"value_attr": "data"},
        "signal_ss_energy": {"value_attr": "num"},
    }

    @property
    def feat(self):
        "Simulate one-to-one relationship between chelatase objects."
        num = self.feat_set.count()
        if num == 0:
            return None
        elif num == 1:
            return self.feat_set.first()
        else:
            raise ValueError("ChelataseFshift '%s' has '%s' feats" % (self, num))

    def create_all_params(self):
        super().create_all_params()
        self._make_poly_a_slippery_site()

    def _make_poly_a_slippery_site(self):
        if self.len != -1:
            return

        if self.feat is None:
            return

        poly_a_site = search_poly_a_slippery(
            self.feat.parent.seq.seq,
            self.feat.parent.strand,
            (self.feat.parent.start, self.feat.parent.end),
            self.len,
            self.org.transl_table,
        )
        if poly_a_site is not None:
            self.set_param("poly_a_start", num=poly_a_site[0])
            self.set_param("poly_a_end", num=poly_a_site[1])

    def _make_signal_secondary_structure(self):
        if self.prm_dict.get("poly_a_start") is None:
            return
        seq = self.seq.seq
        max_length = 200
        if self.strand == 1:
            poly_a_start = self.prm_dict["poly_a_start"]
            seq = seq[poly_a_start : poly_a_start + max_length]
            seq = seq.transcribe()
        elif self.strand == -1:
            poly_a_end = self.prm_dict["poly_a_end"]
            seq = seq[poly_a_end - max_length : poly_a_end]
            seq = seq.reverse_complement().transcribe()
        else:
            raise ValueError
        ss_data = search_stable_secondary_structure(str(seq), 30, max_length)
        rna_seq = ss_data["seq"]
        energy = ss_data["energy"]
        length = len(rna_seq)
        ss_struct = ss_data["pred"].split("\n")[-1]
        if self.strand == 1:
            ss_start = poly_a_start
            ss_end = poly_a_start + length
        elif self.strand == -1:
            ss_end = poly_a_end
            ss_start = poly_a_end - length
        else:
            ValueError
        self.set_param("signal_ss_start", num=ss_start)
        self.set_param("signal_ss_end", num=ss_end)
        self.set_param("signal_ss_struct", data=ss_struct)
        self.set_param("signal_ss_seq", data=rna_seq)
        self.set_param("signal_ss_energy", num=energy)
