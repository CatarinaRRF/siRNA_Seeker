from django.urls import path
from .views import *
from django.contrib.auth import views as auth_views

urlpatterns = [
    # Rota de Cadastro
    path("accounts/signup/", RegisterView.as_view(), name="SignUp"),
    
    # Rota de Login (Utiliza o formulário customizado do Tabler)
    path('accounts/login/', CustomLoginView.as_view(
        redirect_authenticated_user=True, 
        template_name='user/sign_in.html', 
        authentication_form=LoginForm
    ), name='login'),
    
    # Rota de Perfil do Usuário
    path("accounts/profile/", ProfileView.as_view(), name="profile"),
    
    # Rota de Logout
    path('accounts/logout/', auth_views.LogoutView.as_view(), name='logout'),
    
    # NOVA ROTA: Recuperação de Senha Facilitada (Direta no Site)
    # Esta rota agora gerencia tanto a verificação do e-mail quanto a troca da senha
    path('accounts/password-reset/', ResetPasswordView.as_view(), name='password_reset'),
]
