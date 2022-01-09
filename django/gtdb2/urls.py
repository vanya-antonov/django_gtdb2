from django.urls import path
from . import views


urlpatterns = [
    path('', views.OrgListView.as_view(), name='gtdb2-home'),
    path('all_cofs/', views.CofListView.as_view(), name='gtdb2-all-cofs'),
    path('org/<int:pk>/', views.OrgDetailView.as_view(), name='gtdb2-org-detail'),
    path('about/', views.about, name='gtdb2-about'),
]

