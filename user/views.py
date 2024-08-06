from django.shortcuts import render, redirect 
from django.contrib import messages
from django.views import View

from .forms import RegisterForm

# Create your views here.

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

class LogInView(View):

    def get(self, request, *args, **kwargs):
        a = print('insira sua função aqui')
        return render(request, self.template_name, a)
    
class LogOutView(View):
    """ 
    O logout pode ser um botão que redireciona para a 
    pag de login, e no backend sai da conta
    não precisa fazer um template para ele

    """
    def get(self, request, *args, **kwargs):
        a = print('insira sua função aqui')
        return render(request, self.template_name, a)

class ProfileView(View):

    def get(self, request, *args, **kwargs):
        a = print('insira sua função aqui')
        return render(request, self.template_name, a)