from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import ListView, DetailView

from .models import Org, Cof


class OrgListView(ListView):
    model = Org
    # These are Django defaults that can be changed if needed
    # template_name = 'gtdb2/org_list.html'
    # context_object_name = 'object_list'
    ordering = ['name']


class CofListView(ListView):
    model = Cof
    template_name = 'gtdb2/cof_list.html'
    context_object_name = 'object_list'


class OrgDetailView(DetailView):
    model = Org


def about(request):
    return render(request, 'gtdb2/about.html', {'title': 'About'})

