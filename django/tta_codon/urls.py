from django.urls import path
from . import views


urlpatterns = [
    path('', views.TtaOrgListView.as_view(), name='tta-codon-home'),
    path('all_cofs/', views.TtaCofListView.as_view(), name='tta-all-cofs'),
    path('index', views.index, name='index'),
#    path('', views.OrgListView.as_view(), name='gtdb2-home'),
#    path('org/<int:pk>/', views.OrgDetailView.as_view(), name='gtdb2-org-detail'),
#    path('about/', views.about, name='gtdb2-about'),
]

