from functools import lru_cache

from django.db.models import Prefetch, Count, Q
from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse
from django.views.generic import ListView, DetailView
from rest_framework.viewsets import ReadOnlyModelViewSet
from rest_framework.response import Response
from rest_framework_extensions.cache.decorators import (
    cache_response
)
from rest_framework_extensions.cache.mixins import CacheResponseMixin

from chelatase_db.models import ChelataseOrg, ChelataseFshift, ChelataseFeat, ChelataseSeq
from chelatase_db.serializers import (
    ChelataseOrgBaseSerializer,
    ChelataseOrgDetailSerializer,
    ChelataseFshiftSerializer,
    sort_organisms
)
from gtdb2.models import FeatParam, OrgParam


class OrgApiViewSet(ReadOnlyModelViewSet):
    queryset = ChelataseOrg.objects.filter(
        seq__feat__param_set__name="chel_subunit", seq__feat__param_set__value="M"
    ).distinct()

    @cache_response(None)
    def list(self, request):
        queryset = self.get_queryset()
        
        chel_feats_prefetch = ChelataseFeat.objects.prefetch_related(
            Prefetch(
                "param_set",
                queryset=FeatParam.objects.filter(name__in=["chel_gene_group", "chel_evalue"]),
                to_attr="param_prefetch",
            )
        ).filter(param_set__name="chel_gene_group")

        seq_prefetch = ChelataseSeq.objects.prefetch_related(
            Prefetch("feat_set", queryset=chel_feats_prefetch, to_attr="feat_prefetch")
        )
        org_taxonomy_prefetch = OrgParam.objects.filter(name="taxonomy").order_by("num")

        queryset = queryset.prefetch_related(
            Prefetch("seq_set", to_attr="seq_prefetch", queryset=seq_prefetch),
            Prefetch(
            "param_set", to_attr="taxonomy_params", queryset=org_taxonomy_prefetch
        ),
        )
        queryset=queryset.annotate(
            num_fs=Count("seq__feat", filter=Q(seq__feat__type="fsCDS")))
        queryset = sort_organisms(queryset)[::-1]
        serializer = ChelataseOrgBaseSerializer(queryset, many=True)
        return Response(serializer.data)

    def retrieve(self, request, pk=None):
        queryset = self.get_queryset()
        org = get_object_or_404(queryset, pk=pk)
        serializer = ChelataseOrgDetailSerializer(org)
        return Response(serializer.data)


class FshiftsWithSignalStructureViewSet(CacheResponseMixin, ReadOnlyModelViewSet):
    queryset = ChelataseFshift.objects.filter(
        len=-1, feat__param_set__name="chel_subunit", param_set__name="signal_ss_struct",
    ).all()
    serializer_class = ChelataseFshiftSerializer


class OrgListView(ListView):
    model = ChelataseOrg
    # These are Django defaults that can be changed if needed
    template_name = "chelatase_db/org_list.html"
    # context_object_name = 'object_list'
    # ordering = ['name']
    queryset = ChelataseOrg.objects.filter(
        seq__feat__param_set__name="chel_subunit", seq__feat__param_set__value="M"
    ).distinct()


class OrgDetailView(DetailView):
    model = ChelataseOrg
    template_name = "chelatase_db/org_detail.html"

    def get_context_data(self, **kwargs):
        """Extend Django's method to add additional values to context:
        https://docs.djangoproject.com/en/2.2/ref/class-based-views/generic-display/
        """
        context = super().get_context_data(**kwargs)

        # About sorted: https://stackoverflow.com/a/8018989/310453
        b12_dict = self.object.get_full_pathway_gene_dict(pathway="B12")
        chl_dict = self.object.get_full_pathway_gene_dict(pathway="Chlorophyll")
        context["b12_genes"] = sorted(b12_dict.items())
        context["chl_genes"] = sorted(chl_dict.items())

        return context


def about(request):
    return render(request, "chelatase_db/about.html", {"title": "About"})
