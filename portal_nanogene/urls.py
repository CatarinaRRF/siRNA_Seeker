"""
URL configuration for portal_nanogene project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf.urls.static import static
from django.conf import settings 
from django.views.generic import TemplateView

urlpatterns = [
    path("admin/", admin.site.urls),
    path('', include('dashboard.urls')), # Configura o Dashboard
    path('', include('faq.urls')), # Configura paginas de documentação e FAQ
    path('', include('pages.urls')), # Configura paginas staticas do site
    path('', include('sirna_seeker.urls')), # Configura o app siRNA Seeker
    path('', include('user.urls')), # Configura usuarios 

    
    
    path('celery-progress/', include('celery_progress.urls')), #config do celery
]

# Custom /media files  
urlpatterns += staticfiles_urlpatterns()
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

# Custom 404 error view
handler404 = 'pages.views.error_404' 
# Custom 500 error view
handler500 = 'pages.views.error_500' 