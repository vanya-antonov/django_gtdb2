from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import ListView, DetailView

from .models import ChelataseOrg


class OrgListView(ListView):
    model = ChelataseOrg
    # These are Django defaults that can be changed if needed
    template_name = 'chelatase_db/org_list.html'
    # context_object_name = 'object_list'
    ordering = ['name']


class OrgDetailView(DetailView):
    model = ChelataseOrg
    template_name = 'chelatase_db/org_detail.html'


def about(request):
    return render(request, 'chelatase_db/about.html', {'title': 'About'})

