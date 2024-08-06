from django.urls import path
from .views import form_view, loading
from django.views.generic import TemplateView

urlpatterns = [
    path('sirna/', TemplateView.as_view(template_name='sirna/sirna_home.html'), name='sirna_home'), # Pagina principal  
    path('sirna/search/', form_view, name='search'), # Pagina com formulario
    path('sirna/loading/<str:token>/', loading, name='loading'), # Pagina que roda o algoritimo, mostra progresso
]