from django.urls import path
from .views import dashboard, teste

urlpatterns = [
    path('dashboard', dashboard, name='dashboard'),
    path('teste', teste, name='teste')
]