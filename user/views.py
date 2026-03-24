#-------------------------------------------------------------------#
# Django utils
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages
from django.views import View
from django.contrib.auth.models import User
from django.contrib.auth.forms import SetPasswordForm
from django.contrib.auth.views import LoginView
from django.contrib.auth import logout as auth_logout
from django.db import transaction

# Modelos e Formulários Locais
from .models import UserProfile
from .forms import RegisterForm, LoginForm, UserProfileForm
#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
#                        Registo de Utilizadores                    #
#-------------------------------------------------------------------#
class RegisterView(View):
    form_class = RegisterForm
    template_name = 'user/sign_up.html'

    def get(self, request, *args, **kwargs):
        form = self.form_class()
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):
        form = self.form_class(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            messages.success(request, f'Conta criada com sucesso para {username}!')
            return redirect(to='login')
        return render(request, self.template_name, {'form': form})

#-------------------------------------------------------------------#
#                        Login Customizado                          #
#-------------------------------------------------------------------#
class CustomLoginView(LoginView):
    """
    View de Login que suporta "Lembrar-me" e exibe mensagens de erro explícitas.
    """
    form_class = LoginForm

    def form_valid(self, form):
        remember_me = form.cleaned_data.get('remember_me')
        if not remember_me:
            # Se não marcar "lembrar-me", a sessão expira ao fechar o navegador
            self.request.session.set_expiry(0)
            self.request.session.modified = True
        return super(CustomLoginView, self).form_valid(form)

    def form_invalid(self, form):
        """
        Caso o login falhe, adicionamos uma mensagem de erro para garantir
        que o utilizador saiba que a senha ou o usuário estão incorretos.
        """
        messages.error(self.request, "Utilizador ou palavra-passe incorretos. Por favor, tente novamente.")
        return super().form_invalid(form)

#-------------------------------------------------------------------#
#                        Perfil do Utilizador                       #
#-------------------------------------------------------------------#
class ProfileView(View):
    template_name = 'user/profile.html'

    def get(self, request):
        from django_celery_results.models import TaskResult
        user_profile = get_object_or_404(UserProfile, user=request.user)
        form = UserProfileForm(instance=user_profile)
        tasks = TaskResult.objects.all()
        return render(request, self.template_name, {'form': form, 'tasks': tasks, 'user_profile': user_profile})

    def post(self, request):
        from django_celery_results.models import TaskResult
        user_profile = get_object_or_404(UserProfile, user=request.user)
        form = UserProfileForm(request.POST, instance=user_profile)
        
        # Lógica para apagar tarefas do Celery
        if 'delete_task' in request.POST:
            task_id = request.POST.get('task_id')
            task = get_object_or_404(TaskResult, task_id=task_id)
            task.delete()
            messages.success(request, "Tarefa eliminada com sucesso.")
            return redirect('profile')

        if form.is_valid():
            form.save()
            messages.success(request, "Perfil atualizado!")
            return redirect('profile')

        return render(request, self.template_name, {'form': form, 'user_profile': user_profile})

#-------------------------------------------------------------------#
#             Recuperação de Palavra-passe Ultra-Segura             #
#-------------------------------------------------------------------#
class ResetPasswordView(View):
    """
    View para redefinição de senha sem dependência de e-mail.
    Garante a invalidação de sessões antigas e limpeza de cache de CSRF.
    Lida com múltiplos usuários que possuem o mesmo e-mail.
    """
    template_name = 'user/forget_password.html'

    def get(self, request):
        # Limpeza total ao carregar a página para evitar conflitos de CSRF
        auth_logout(request)
        request.session.flush()
        return render(request, self.template_name, {'step': 1})

    def post(self, request):
        # PASSO 1: Verificar se o e-mail existe no sistema
        if 'check_email' in request.POST:
            email = request.POST.get('email', '').strip()
            # Busca ignorando maiúsculas/minúsculas
            users = User.objects.filter(email__iexact=email)
            
            if users.exists():
                # Se houver duplicados, pegamos o primeiro, mas avisamos o username
                user = users.first()
                request.session['reset_user_id'] = user.id
                form = SetPasswordForm(user)
                
                # Mensagem importante para o usuário saber QUAL conta está mexendo
                msg = f"Iniciando recuperação para o usuário: '{user.username}'."
                if users.count() > 1:
                    usernames = ", ".join([u.username for u in users])
                    msg = f"Atenção: Existem {users.count()} contas com este e-mail ({usernames}). Recuperando acesso apenas para: '{user.username}'."
                
                messages.warning(request, msg)
                return render(request, self.template_name, {
                    'form': form, 
                    'step': 2, 
                    'target_user': user.username
                })
            else:
                messages.error(request, f"O e-mail '{email}' não foi encontrado.")
                return render(request, self.template_name, {'step': 1})

        # PASSO 2: Validar e salvar a nova palavra-passe
        if 'set_password' in request.POST:
            user_id = request.session.get('reset_user_id')
            if not user_id:
                messages.error(request, "Sessão expirada. Por favor, recomece.")
                return redirect('password_reset')
            
            user = get_object_or_404(User, id=user_id)
            form = SetPasswordForm(user, request.POST)
            
            if form.is_valid():
                # Usamos uma transação para garantir que a gravação é atómica
                with transaction.atomic():
                    # Definimos a senha manualmente para garantir que o set_password é chamado
                    new_password = form.cleaned_data.get('new_password1')
                    user.set_password(new_password)
                    user.save()
                
                # Mensagem de sucesso com o Username para o usuário não esquecer
                messages.success(request, f"Sucesso! Palavra-passe alterada para o utilizador: {user.username}")
                
                # LIMPEZA TOTAL DA SESSÃO
                # Crucial para o botão de login funcionar à primeira tentativa
                auth_logout(request)
                request.session.flush()
                
                return redirect('login')
            
            # Se houver erros de validação no formulário
            return render(request, self.template_name, {
                'form': form, 
                'step': 2, 
                'target_user': user.username
            })

        return redirect('password_reset')
