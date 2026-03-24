#-------------------------------------------------------------------#
# Django utils
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages
from django.views import View
from django.contrib.auth.models import User
from django.contrib.auth.forms import SetPasswordForm
from django.contrib.auth.views import LoginView
from django.contrib.auth import update_session_auth_hash

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
    form_class = LoginForm

    def form_valid(self, form):
        remember_me = form.cleaned_data.get('remember_me')
        if not remember_me:
            # Se não marcar "lembrar-me", a sessão expira ao fechar o navegador
            self.request.session.set_expiry(0)
            self.request.session.modified = True
        return super(CustomLoginView, self).form_valid(form)

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
#             Recuperação de Palavra-passe Corrigida                #
#-------------------------------------------------------------------#
class ResetPasswordView(View):
    """
    View para redefinição de senha sem dependência de envio de e-mail.
    Garante a invalidação de sessões antigas e limpeza de cache de CSRF.
    """
    template_name = 'user/forget_password.html'

    def get(self, request):
        # Limpa qualquer resquício de tentativa anterior ao carregar a página
        if 'reset_user_id' in request.session:
            del request.session['reset_user_id']
        return render(request, self.template_name, {'step': 1})

    def post(self, request):
        # PASSO 1: Verificar se o e-mail existe no sistema
        if 'check_email' in request.POST:
            email = request.POST.get('email')
            users = User.objects.filter(email=email)
            
            if users.exists():
                user = users.first()
                # Armazena o ID do utilizador na sessão de forma temporária
                request.session['reset_user_id'] = user.id
                form = SetPasswordForm(user)
                return render(request, self.template_name, {
                    'form': form, 
                    'step': 2, 
                    'target_user': user.username
                })
            else:
                messages.error(request, "E-mail não encontrado no sistema.")
                return render(request, self.template_name, {'step': 1})

        # PASSO 2: Validar e salvar a nova palavra-passe
        if 'set_password' in request.POST:
            user_id = request.session.get('reset_user_id')
            if not user_id:
                messages.error(request, "Sessão expirada. Por favor, recomece o processo.")
                return redirect('password_reset')
            
            user = get_object_or_404(User, id=user_id)
            form = SetPasswordForm(user, request.POST)
            
            if form.is_valid():
                # 1. Salva a nova senha (faz o hash internamente)
                user = form.save()
                
                # 2. Invalida todas as sessões antigas
                # Isto impede que a senha antiga continue a funcionar em outros dispositivos
                update_session_auth_hash(request, user)
                
                # 3. Mensagem informativa com o Username correto
                messages.success(request, f"Sucesso! Palavra-passe alterada para o utilizador: {user.username}")
                
                # 4. LIMPEZA TOTAL DA SESSÃO
                # Resolve o problema do clique duplo no login ao forçar um novo token CSRF
                request.session.flush() 
                
                return redirect('login')
            
            # Se o formulário tiver erros (ex: senhas não coincidem)
            return render(request, self.template_name, {
                'form': form, 
                'step': 2, 
                'target_user': user.username
            })

        return redirect('password_reset')
