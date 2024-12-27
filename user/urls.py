from django.urls import path
from .views import *
from django.contrib.auth import views as auth_views
from django.views.generic import TemplateView


urlpatterns = [
    path("accounts/signup/", RegisterView.as_view(), name="SignUp"), # formata o singup
    path('accounts/login/', CustomLoginView.as_view(redirect_authenticated_user=True, template_name='user/sign_in.html', authentication_form=LoginForm), name='login'),
    path("accounts/profile/", ProfileView.as_view(), name="profile"), #formata o profile
    path('accounts/logout/', auth_views.LogoutView.as_view(), name='logout'),
    path('accounts/password-reset/', ResetPasswordView.as_view(), name='password_reset'), #reset
    path('accounts/password-reset-confirm/<uidb64>/<token>/',
         auth_views.PasswordResetConfirmView.as_view(template_name='user/password_reset_confirm.html'),
         name='password_reset_confirm'), #confirm
    path(
        'accounts/password-reset-complete/',
        auth_views.PasswordResetCompleteView.as_view(template_name='user/password_reset_complete.html'),
        name='password_reset_complete'), #complete
]

