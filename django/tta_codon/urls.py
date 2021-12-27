from django.urls import path
from . import views


urlpatterns = [
    path('', views.TtaOrgListView.as_view(), name='tta-codon-home'),
    path('all_cofs/', views.TtaCofListView.as_view(), name='tta-all-cofs'),
    path('cof/<int:pk>/', views.TtaCofDetailView.as_view(), name='tta-cof-detail'),
    path('index', views.index, name='index'),
#    path('', views.OrgListView.as_view(), name='gtdb2-home'),
#    path('about/', views.about, name='gtdb2-about'),
]

