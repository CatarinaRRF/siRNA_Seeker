# Libraries
import uuid
from django.urls import reverse
from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from django.http import Http404
import os
import logging

# Models and forms
from .forms import form_search
from .models import model_forms
from django_celery_results.models import TaskResult

# Directories
from .algorithm import *
from .tasks import *

from user.models import UserProfile

# Debugging the functions
logger = logging.getLogger(__name__)


#-------------------------------------------------------------------------------#
#                                   Forms View
#-------------------------------------------------------------------------------#
@login_required(redirect_field_name='next', login_url="/accounts/login/")
def form_view(request):
    user = request.user
    user_profile, _ = UserProfile.objects.get_or_create(user=user)


    if request.method == "POST":
        form = form_search(request.POST, request.FILES)
        if form.is_valid():
            form_instance = form.save(commit=False)
            form_instance.user = request.user
            form_instance.save()
            
            # Debugging
            print("Run Blast:", form.cleaned_data['run_blast'])  
            logger.info("Run Blast: %s", form.cleaned_data['run_blast'])
            
            # Gerar um token único
            token = str(uuid.uuid4())
            
            # Armazenar o token na sessão
            request.session['loading_token'] = token
            
            # Redirecionar para a página 'loading' com o token
            return redirect(reverse('loading', kwargs={'token': token}))
        
        else:
            logger.error("Form is not valid: %s", form.errors)
    else:
        form = form_search()
        # Adiciona informações do usuário ao contexto
    context = {
        'form': form,
        'user_profile': {
            'name': user.get_full_name(),
            'profission': user_profile.profission,
            'profile_image': user_profile.profile_image,
        }
    }
    return render(request, 'sirna/forms.html', context)


#-------------------------------------------------------------------------------#
#                                Run the alghorithm
#-------------------------------------------------------------------------------#
@login_required(redirect_field_name='next', login_url="/accounts/login/")
def loading(request, token):
    user_p = request.user
    user_profile, _ = UserProfile.objects.get_or_create(user=user_p)
    # Obtém o token da sessão
    session_token = request.session.get('loading_token')

    # Verifica se o token é válido
    if session_token != token:
        # Se o token não corresponder, redireciona para uma página de erro ou inicial
        raise Http404("Not Found")

    # Remove o token da sessão após verificação
    del request.session['loading_token']

    model_forms_data = model_forms.objects.filter(user=request.user).last()
    
    _directory = 'media'

    # User
    user = request.user.username

    # Criando o path para a sequência
    sequence = str(model_forms_data.sequence)
    sequence_path = os.path.join(_directory, sequence)

    # Armazenando outras variáveis
    sequence_tag = str(model_forms_data.sequence_tag)
    autor = model_forms_data.autor
    size = int(model_forms_data.size)
    include_tm = bool(model_forms_data.include_tm)
    max_tm = float(model_forms_data.max_tm)
    run_blast = bool(model_forms_data.run_blast)
    organism = model_forms_data.organism
    database = model_forms_data.database
    identity = model_forms_data.identity
    query_cover = model_forms_data.query_cover

    # Executar tarefa Celery se houver um registro no banco de dados
    if model_forms_data:
        task = selection.delay(user, sequence_path, sequence_tag, autor,
                               size, include_tm, max_tm, run_blast, organism,
                               database, identity, query_cover)
        try:
            task_result = TaskResult.objects.get(task_id=task.task_id)
            task_result.task_name = sequence_tag  # Atualizando o nome da tarefa
            task_result.save()
        except TaskResult.DoesNotExist:
            pass

    context = {
        'task_id': task.task_id,
        'user_profile': {
            'name': user_p.get_full_name(),
            'profission': user_profile.profission,
            'profile_image': user_profile.profile_image,
        }
    }

    return render(request, 'sirna/loadings.html', context)

            
