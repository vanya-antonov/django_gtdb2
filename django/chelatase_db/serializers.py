from rest_framework import serializers
from chelatase_db.models.org import ChelataseOrg, ChelataseFeat


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


class ChelataseOrgBaseSerializer(serializers.ModelSerializer):
    genotype = serializers.CharField(source="prm_dict.chel_genotype_genes")

    class Meta:
        model = ChelataseOrg
        fields = ("id", "name", "phylum", "kingdom", "genotype")


class ChelataseOrgDetailSerializer(serializers.ModelSerializer):
    chel_subunit_feat_set = serializers.ListField(child=ChelataseFeatBaseSerializer())
    b12_genes_count = serializers.SerializerMethodField()
    chl_genes_count = serializers.SerializerMethodField()

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
        )

    def get_b12_genes_count(self, obj):
        b12_dict = obj.get_full_pathway_gene_dict(pathway="B12")
        return sorted([(k,len(v))for k,v in b12_dict.items()])

    def get_chl_genes_count(self, obj):
        chl_dict = obj.get_full_pathway_gene_dict(pathway="Chlorophyll")
        return sorted([(k,len(v))for k,v in chl_dict.items()])
