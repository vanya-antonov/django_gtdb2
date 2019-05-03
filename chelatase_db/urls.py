from django.urls import path
from . import views


urlpatterns = [
    path('', views.OrgListView.as_view(), name='chelatase-home'),
    path('org/<int:pk>/', views.OrgDetailView.as_view(),
         name='chelatase-org-detail'),
    path('about/', views.about, name='chelatase-about'),
]

