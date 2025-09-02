from django.urls import path
from django.views.generic import TemplateView
from . import views

urlpatterns = [
    # Static Pages
    path('', TemplateView.as_view(template_name='pages/main.html'), name='NanoGen'), #Pagina principal
    #path('projects', views.project_list, name='project_list'), #Projects 
    path('projects', TemplateView.as_view(template_name='pages/projects_test.html'), name='project_list'), #Projects 
    path('authors', TemplateView.as_view(template_name='pages/god.html'), name='authors'), #Events
    path('events', TemplateView.as_view(template_name='pages/events.html'), name='events'), #Events
    path('tools', TemplateView.as_view(template_name='pages/tools.html'), name='Tools'), #Pagina que coleta todos as ferramentas do grupo # Pagina com as informações de contato 
    
    #Reaserch Groups
    path('armando', TemplateView.as_view(template_name='clowns/armando.html'), name='armando'), 
    path('frederico', TemplateView.as_view(template_name='clowns/fredd.html'), name='frederico'), 
    path('luciana', TemplateView.as_view(template_name='clowns/luciana.html'), name='luciana'), 
    path('mauro', TemplateView.as_view(template_name='clowns/mauro.html'), name='mauro'), 
    path('silvia', TemplateView.as_view(template_name='clowns/silvia.html'), name='silvia'), 
    path('matheus', TemplateView.as_view(template_name='clowns/matheus.html'), name='matheus'), 

    
    path('erro_404', TemplateView.as_view(template_name='errors/error_404.html'), name='erro_404'),
    #path('authors', views.authors, name='authors'),
]