from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse
from django.views.generic import ListView, DetailView
from rest_framework.viewsets import ReadOnlyModelViewSet
from rest_framework.response import Response

from chelatase_db.models import ChelataseOrg, ChelataseFshift
from chelatase_db.serializers import (
    ChelataseOrgBaseSerializer,
    ChelataseOrgDetailSerializer,
    ChelataseFshiftSerializer,
)


class OrgApiViewSet(ReadOnlyModelViewSet):
    queryset = ChelataseOrg.objects.filter(
        seq__feat__param_set__name="chel_subunit", seq__feat__param_set__value="M"
    ).distinct()

    def list(self, request):
        queryset = self.get_queryset()
        serializer = ChelataseOrgBaseSerializer(queryset, many=True)
        return Response(serializer.data)

    def retrieve(self, request, pk=None):
        queryset = self.get_queryset()
        org = get_object_or_404(queryset, pk=pk)
        serializer = ChelataseOrgDetailSerializer(org)
        return Response(serializer.data)


class FshiftsWithSignalStructureViewSet(ReadOnlyModelViewSet):
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
