from Bio.Seq import Seq
from django.test import TestCase

from chelatase_db.models.fshift import search_poly_a_slippery


class TestFindSlippery(TestCase):
    def test_search_poly_a_slippery(self):
        transl_table = 11  # The Bacterial, Archaeal and Plant Plastid Code
        fshift_len = -1  # poly-a for -1 only!
        seq_str = "cgtgctcgatatgtcgtgctcgatcatctaccggattaaggcaccggactctagcccaggggcaaaaaagccggtcaaaagtaggattgttgcaaaaacctcagcttaggcaaaagccg"
        dna = Seq(seq_str)
        span_1 = search_poly_a_slippery(dna, 1, (10, 109), fshift_len, transl_table)
        print(span_1)
        self.assertEqual(span_1, (93, 98))
        span_2 = search_poly_a_slippery(dna, 1, (10, 84), fshift_len, transl_table)
        self.assertEqual(span_2, (63, 69))
        span_3 = search_poly_a_slippery(dna, 1, (66, 84), fshift_len, transl_table)
        print(span_3)
        self.assertIsNone(span_3, None)
        span_1_reverce = search_poly_a_slippery(
            dna.reverse_complement(), -1, (len(dna) - 109, len(dna) - 10), fshift_len, transl_table,
        )
        self.assertEqual(span_1_reverce, (len(dna) - span_1[1], len(dna) - span_1[0]))
        span_2_reverce = search_poly_a_slippery(
            dna.reverse_complement(), -1, (len(dna) - 84, len(dna) - 10), fshift_len, transl_table,
        )
        self.assertEqual(span_2_reverce, (len(dna) - span_2[1], len(dna) - span_2[0]))
        span_3_reverce = search_poly_a_slippery(
            dna.reverse_complement(), -1, (len(dna) - 84, len(dna) - 66), fshift_len, transl_table,
        )
        self.assertIsNone(span_3_reverce, None)
