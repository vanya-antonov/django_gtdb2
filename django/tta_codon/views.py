from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import ListView, DetailView

from .models import TtaOrg


class TtaOrgListView(ListView):
    model = TtaOrg
    template_name = "tta_codon/org_list.html"
    # These are Django defaults that can be changed if needed
    # template_name = 'gtdb2/org_list.html'
    # context_object_name = 'object_list'
    ordering = ['name']
    queryset = TtaOrg.objects.filter(name__startswith='Streptomyces')


# Create your views here.
def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")


