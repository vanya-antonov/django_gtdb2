from collections import defaultdict

import math
from rest_framework import serializers
from chelatase_db.models.org import ChelataseOrg, ChelataseFeat, ChelataseFshift

ORGANISMS_HEATMAP_TAXAS = {
    "Proteobacteria": 0,
    "Actinobacteria": 1,
    "Archaea": 2,
    "Chloroflexi": 3,
    "Firmicutes": 4,
    "Cyanobacteria": 5,
    "Other": 6,
}


def taxonomy_to_heatmap_taxa(taxonomy):
    #     print(taxonomy)
    if taxonomy[0] == "Archaea":
        return "Archaea"
    if taxonomy[0] == "Bacteria":
        bacteria_taxa = taxonomy[1]
        if bacteria_taxa in ORGANISMS_HEATMAP_TAXAS:
            return bacteria_taxa
    return "Other"


def log_evalues(evalue):
    return 100 if evalue < 1e-100 else -math.log10(evalue)


def sort_organisms(organisms):

    organisms_taxonomy = []
    for org in organisms:
        taxonomy = tuple(param.value for param in org.taxonomy_params)
        heatmap_taxa_num = ORGANISMS_HEATMAP_TAXAS[taxonomy_to_heatmap_taxa(taxonomy)]
        organisms_taxonomy.append((org, (heatmap_taxa_num, taxonomy)))
    organisms_sorted, _ = zip(*sorted(organisms_taxonomy, key=lambda a: a[1]))
    return organisms_sorted


class ChelataseFeatBaseSerializer(serializers.ModelSerializer):
    frameshift_len = serializers.IntegerField(source="fshift.len", allow_null=True)

    class Meta:
        model = ChelataseFeat
        fields = ("name", "frameshift_len", "prm_dict", "seq_id", "start", "end", "strand", "descr")

    def to_representation(self, instance):
        ret = super().to_representation(instance)
        prm_dict = ret["prm_dict"]
        del prm_dict["seq_nt"]
        translation = prm_dict.pop("translation")
        prm_dict["translation_len"] = len(translation)
        ret["prm_dict"] = prm_dict
        return ret


class ChelataseFshiftSerializer(serializers.ModelSerializer):
    org_name = serializers.CharField(source="org.name")
    org_id = serializers.IntegerField(source="org.id")
    poly_a_slippery = serializers.SerializerMethodField()
    signal = serializers.SerializerMethodField()

    class Meta:
        model = ChelataseFshift
        fields = ("org_name", "org_id", "name", "poly_a_slippery", "signal", "strand")

    def get_poly_a_slippery(self, obj):
        prm_dict = obj.prm_dict
        if "poly_a_start" not in prm_dict:
            return None
        return (prm_dict["poly_a_start"], prm_dict["poly_a_end"])

    def get_signal(self, obj):
        prm_dict = obj.prm_dict
        if "signal_ss_struct" not in prm_dict:
            return None
        data = {}
        data["struct"] = prm_dict["signal_ss_struct"]
        data["seq"] = prm_dict["signal_ss_seq"]
        if obj.strand == 1:
            poly_a_coord = (
                -prm_dict["signal_ss_start"] + prm_dict["poly_a_start"],
                -prm_dict["signal_ss_start"] + prm_dict["poly_a_end"],
            )
            parent_feat_end = obj.feat.parent.end - prm_dict["signal_ss_start"]
            stop_codon_coord = (parent_feat_end - 3, parent_feat_end)
        if obj.strand == -1:
            poly_a_coord = (
                prm_dict["signal_ss_end"] - prm_dict["poly_a_end"],
                prm_dict["signal_ss_end"] - prm_dict["poly_a_start"],
            )
            parent_feat_end = prm_dict["signal_ss_end"] - obj.feat.parent.start
            stop_codon_coord = (parent_feat_end - 3, parent_feat_end)
        data["poly_a_coord"] = poly_a_coord
        data["stop_codon_coord"] = stop_codon_coord
        data["energy"] = prm_dict["signal_ss_energy"]
        return data


class ChelataseOrgBaseSerializer(serializers.ModelSerializer):
    genotype = serializers.SerializerMethodField()
    evalues = serializers.SerializerMethodField()
    heatmap_taxa = serializers.SerializerMethodField()
    num_fshifts = serializers.IntegerField(source='num_fs')

    gene_groups = (
        "chlI_bchI",
        "chlD_bchD",
        "chlH_bchH",
        "cobN",
        "cobS",
        "cobT",
        "bchE",
        "chlB_bchB",
        "chlG_bchG",
        "chlL_bchL",
        "chlM_bchM",
        "chlN_bchN",
        "cobD_cobC",
        "cobO",
        "cobP_cobU",
        "cobQ",
        "cobV_cobS",
        "cysG_cobA",
    )

    class Meta:
        model = ChelataseOrg
        fields = (
            "id",
            "name",
            "phylum",
            "kingdom",
            "genotype",
            "evalues",
            "heatmap_taxa",
            "num_fshifts",
        )

    def get_genotype(self, obj):
        return obj.prm_dict['chel_genotype_genes'].split(', ')

    def get_evalues(self, obj):
        chel_evalues = defaultdict(list)
        for seq in obj.seq_prefetch:
            for feat in seq.feat_prefetch:
                if feat.pare
                chelevals = {param.name: param.value for param in feat.param_prefetch}
                chel_evalues[chelevals["chel_gene_group"]].append(float(chelevals["chel_evalue"]))
        chel_evalues = {k: min(v) for k, v in chel_evalues.items() if k in self.gene_groups}
        gene_groupchel_evalues = {
            gene_group: log_evalues(chel_evalues.get(gene_group, 1))
            for gene_group in self.gene_groups
        }
        return gene_groupchel_evalues

    def get_heatmap_taxa(self, obj):
        taxonomy = tuple(param.value for param in obj.taxonomy_params)
        return taxonomy_to_heatmap_taxa(taxonomy)


class ChelataseOrgDetailSerializer(serializers.ModelSerializer):
    chel_subunit_feat_set = serializers.ListField(child=ChelataseFeatBaseSerializer())
    b12_genes_count = serializers.SerializerMethodField()
    chl_genes_count = serializers.SerializerMethodField()
    fshift_set = serializers.ListField(child=ChelataseFshiftSerializer())

    class Meta:
        model = ChelataseOrg
        fields = (
            "id",
            "name",
            "phylum",
            "kingdom",
            "prm_dict",
            "chel_subunit_feat_set",
            "b12_genes_count",
            "chl_genes_count",
            "fshift_set",
        )

    def get_b12_genes_count(self, obj):
        b12_dict = obj.get_full_pathway_gene_dict(pathway="B12")
        return sorted([(k, len(v)) for k, v in b12_dict.items()])

    def get_chl_genes_count(self, obj):
        chl_dict = obj.get_full_pathway_gene_dict(pathway="Chlorophyll")
        return sorted([(k, len(v)) for k, v in chl_dict.items()])
