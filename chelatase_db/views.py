from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import ListView, DetailView

from chelatase_db.models import ChelataseOrg


class OrgListView(ListView):
    model = ChelataseOrg
    # These are Django defaults that can be changed if needed
    template_name = 'chelatase_db/org_list.html'
    # context_object_name = 'object_list'
    # ordering = ['name']
    queryset = ChelataseOrg.objects.filter(
        param_set__name='chel_num_chld',
        param_set__num__gt=0)


class OrgDetailView(DetailView):
    model = ChelataseOrg
    template_name = 'chelatase_db/org_detail.html'

    def get_context_data(self, **kwargs):
        """Extend Django's method to add additional values to context:
        https://docs.djangoproject.com/en/2.2/ref/class-based-views/generic-display/
        """
        context = super().get_context_data(**kwargs)

        # About sorted: https://stackoverflow.com/a/8018989/310453
        b12_dict = self.object.get_full_pathway_gene_dict(pathway='B12')
        chl_dict = self.object.get_full_pathway_gene_dict(pathway='Chlorophyll')
        context['b12_genes'] = sorted(b12_dict.items())
        context['chl_genes'] = sorted(chl_dict.items())

        return context


def about(request):
    return render(request, 'chelatase_db/about.html', {'title': 'About'})

