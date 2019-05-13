from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import ListView, DetailView

from chelatase_db.models import ChelataseOrg
from chelatase_db.lib.config import read_pathway_gene_info


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

        # Get statistics about B12/CHL genes in the given org
        info_dict = read_pathway_gene_info()
        all_gene_groups = {}
        for gene in info_dict.values():
            # e.g. 'cobD_cobC'
            gene_group = gene['chel_gene_group']
            if gene_group not in all_gene_groups:
                all_gene_groups[gene_group] = gene

        context['b12_groups'] = []
        context['chlorophyll_groups'] = []
        for group in all_gene_groups.values():
            if group['chel_pathway'] == 'B12':
                context['b12_groups'].append(group)
            elif group['chel_pathway'] == 'Chlorophyll':
                context['chlorophyll_groups'].append(group)

        return context


def about(request):
    return render(request, 'chelatase_db/about.html', {'title': 'About'})

