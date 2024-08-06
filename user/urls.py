from django.urls import path
from django.views.generic import TemplateView
from .views import *

urlpatterns = [
    path("accounts/signup/", RegisterView.as_view(), name="SignUp"), # formata o singup
    path("accounts/login/", RegisterView.as_view(), name="login"), #formata o login
    path("accounts/profile/", RegisterView.as_view(), name="profile"), #formata o profile
    path("accounts/logout/", RegisterView.as_view(), name="logout"), #formata o logout
]