from django.urls import path
from django.views.generic import TemplateView
from .views import *

urlpatterns = [
    path("accounts/signup/", RegisterView.as_view(), name="SignUp"), # formata o singup
    path("accounts/login/", LogInView.as_view(), name="login"), #formata o login
    path("accounts/profile/", ProfileView.as_view(), name="profile"), #formata o profile
]

