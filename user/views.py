#-------------------------------------------------------------------#
# Django utils
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages
from django.views import View
from django.contrib.auth.models import User
from django.contrib.auth.forms import SetPasswordForm
from django.contrib.auth.views import LoginView
from django.contrib.auth import logout as auth_logout

# Modelos e Formulários Locais
from .models import UserProfile
from .forms import RegisterForm, LoginForm, UserProfileForm
#-------------------------------------------------------------------#

# ... (RegisterView e ProfileView permanecem as mesmas) ...

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

class CustomLoginView(LoginView):
    form_class = LoginForm
    def form_valid(self, form):
        remember_me = form.cleaned_data.get('remember_me')
        if not remember_me:
            self.request.session.set_expiry(0)
            self.request.session.modified = True
        return super(CustomLoginView, self).form_valid(form)

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
#             Recuperação de Palavra-passe Blindada                #
#-------------------------------------------------------------------#
class ResetPasswordView(View):
    template_name = 'user/forget_password.html'

    def get(self, request):
        # Limpa tudo antes de começar para evitar conflitos de CSRF
        auth_logout(request)
        request.session.flush()
        return render(request, self.template_name, {'step': 1})

    def post(self, request):
        # PASSO 1: Verificar e-mail
        if 'check_email' in request.POST:
            email = request.POST.get('email')
            users = User.objects.filter(email=email)
            
            if users.exists():
                user = users.first()
                request.session['reset_user_id'] = user.id
                form = SetPasswordForm(user)
                return render(request, self.template_name, {
                    'form': form, 
                    'step': 2, 
                    'target_user': user.username
                })
            else:
                messages.error(request, "E-mail não encontrado.")
                return render(request, self.template_name, {'step': 1})

        # PASSO 2: Salvar nova senha
        if 'set_password' in request.POST:
            user_id = request.session.get('reset_user_id')
            if not user_id:
                return redirect('password_reset')
            
            user = get_object_or_404(User, id=user_id)
            form = SetPasswordForm(user, request.POST)
            
            if form.is_valid():
                # Forçamos a gravação manual para garantir que o SQLite receba o dado
                new_password = form.cleaned_data.get('new_password1')
                user.set_password(new_password)
                user.save()
                
                # Mensagem de sucesso clara
                messages.success(request, f"Senha atualizada para o usuário: {user.username}. Faça login agora.")
                
                # LIMPAMOS A SESSÃO TOTALMENTE PARA O LOGIN FUNCIONAR NO PRIMEIRO CLIQUE
                auth_logout(request)
                request.session.flush()
                
                return redirect('login')
            
            return render(request, self.template_name, {
                'form': form, 
                'step': 2, 
                'target_user': user.username
            })

        return redirect('password_reset')
