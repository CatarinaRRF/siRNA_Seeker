from django.urls import path
from django.views.generic import TemplateView
from . import views

urlpatterns = [
    # Static Pages
    path('', TemplateView.as_view(template_name='pages/main.html'), name='NanoGen'), #Pagina principal
    path('projects', views.project_list, name='project_list'), #Projects 
    path('events', TemplateView.as_view(template_name='pages/events.html'), name='events'), #Events
    path('tools', TemplateView.as_view(template_name='pages/tools.html'), name='Tools'), #Pagina que coleta todos as ferramentas do grupo # Pagina com as informações de contato

  
]