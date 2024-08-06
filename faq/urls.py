from django.urls import path
from django.views.generic import TemplateView

urlpatterns = [

    path('faq', TemplateView.as_view(template_name='faq/faq.html'), name='faq'), # Perguntas bobas
    path('documentation', TemplateView.as_view(template_name='faq/documentation.html'), name='documentation'), # Explicar como funciona as tools
]