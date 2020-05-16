from django.urls import path, include
from rest_framework.routers import DefaultRouter

from . import views


router = DefaultRouter()
router.register(r'organisms', views.OrgApiViewSet,  basename='org')
api_patterns = router.urls


# api_patterns = [
#     path("orgs", views.OrgListApiView.as_view()),
# ]


urlpatterns = [
    path('', views.OrgListView.as_view(), name='chelatase-home'),
    path('org/<int:pk>/', views.OrgDetailView.as_view(),
         name='chelatase-org-detail'),
    path('about/', views.about, name='chelatase-about'),
    path('api/', include(api_patterns))
]
