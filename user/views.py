#-------------------------------------------------------------------#
# Django utils
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages
from django.views import View
from django.http import JsonResponse
from django.utils.decorators import method_decorator
from django.urls import reverse_lazy

#Django auth
from django.contrib.auth.decorators import login_required
from django.contrib.auth.views import LoginView, PasswordResetView
from django.contrib.messages.views import SuccessMessageMixin
import json

#models and forms
from .models import UserProfile
from .forms import RegisterForm, LoginForm, UserProfileForm
from celery.result import AsyncResult
from django_celery_results.models import TaskResult 
#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
#                        Registration                               #
#-------------------------------------------------------------------#
class RegisterView(View):
    form_class = RegisterForm
    initial = {'key': 'value'}
    template_name = 'user/sign_up.html'

    def get(self, request, *args, **kwargs):
        form = self.form_class(initial=self.initial)
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):
        form = self.form_class(request.POST)

        if form.is_valid():
            form.save()

            username = form.cleaned_data.get('username')
            messages.success(request, f'Account created for {username}')

            return redirect(to='/')

        return render(request, self.template_name, {'form': form})

#-------------------------------------------------------------------#
#                           Login                                   #
#-------------------------------------------------------------------#
class CustomLoginView(LoginView):
    form_class = LoginForm

    def form_valid(self, form):
        remember_me = form.cleaned_data.get('remember_me')

        if not remember_me:
            # set session expiry to 0 seconds. So it will automatically close the session after the browser is closed.
            self.request.session.set_expiry(0)

            # Set session as modified to force data updates/cookie to be saved.
            self.request.session.modified = True

        # else browser session will be as long as the session cookie time "SESSION_COOKIE_AGE" defined in settings.py
        return super(CustomLoginView, self).form_valid(form)

#-------------------------------------------------------------------#
#                          Profile                                  #
#-------------------------------------------------------------------#

class ProfileView(View):
    @method_decorator(login_required)
    def get(self, request, *args, **kwargs):
        user = request.user
        user_profile, _ = UserProfile.objects.get_or_create(user=user)
        form = UserProfileForm(instance=user.profile)

        # Lista de IDs de imagens de perfil
        avatar_ids = list(range(1, 17))

        # Fetch tasks related to the user
        user_tasks = TaskResult.objects.filter(task_creator=user).order_by('-date_created')

        tasks = []
        for task in user_tasks:
            task_id = task.task_id
            task_status = AsyncResult(task_id).status
            try:
                task_result_data = json.loads(task.result)
                query_meta_data = {}
                if isinstance(task_result_data, list) and len(task_result_data) > 1 and isinstance(task_result_data[1], dict):
                    query_meta_data = task_result_data[1]
                task_name = query_meta_data.get('query_title', '') or f"Unnamed"
            except (json.JSONDecodeError, IndexError, TypeError):
                task_name = f"Unnamed"
            tasks.append({
                'id': task.id,
                'task_id': task_id,
                'name': task_name,
                'status': task_status,
                'result': task.result,
                'created': task.date_created,
                'is_running': task_status in ['PENDING', 'STARTED', 'PROGRESS'],
            })

        context = {
            'tasks': tasks,
            'form': form,
            'user_profile': {
                'name': user.get_full_name(),
                'email': user.email,
                'account_anni': user.date_joined,
                'profission': user_profile.profission,
                'institution': user_profile.institution,
                'profile_image': user_profile.profile_image,
            },
        }
        return render(request, 'user/profile.html', context)

    @method_decorator(login_required)
    def post(self, request, *args, **kwargs):
        user = request.user
        user_profile, _ = UserProfile.objects.get_or_create(user=user)
        form = UserProfileForm(request.POST, instance=user_profile)
            
        # Verifica se é para deletar uma tarefa
        if 'delete_task' in request.POST:
            task_id = request.POST.get('task_id')
            task = get_object_or_404(TaskResult, task_id=task_id)
            task.delete()
            messages.success(request, "Task deleted successfully.")
            return redirect('profile')

        # Atualiza o perfil do usuário
        profile_image_path = request.POST.get('profile_image')
        if profile_image_path:
            user_profile.profile_image = profile_image_path

        if form.is_valid():
            form.save()
            user_profile.save()
            return redirect('profile')

        return JsonResponse({'success': False, 'errors': form.errors})

#-------------------------------------------------------------------#
#                        Change Password                            #
#-------------------------------------------------------------------#
class ResetPasswordView(SuccessMessageMixin, PasswordResetView):
    template_name = 'user/forget_password.html'
    email_template_name = 'user/password_reset_email.html'
    subject_template_name = 'user/password_reset_subject.html'
    success_message = "We've emailed you instructions for setting your password, " \
                      "if an account exists with the email you entered. You should receive them shortly." \
                      " If you don't receive an email, " \
                      "please make sure you've entered the address you registered with, and check your spam folder."
    success_url = reverse_lazy('NanoGen')


