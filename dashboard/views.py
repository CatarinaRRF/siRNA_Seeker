
# Django lovries
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.contrib.auth.decorators import login_required
from django_celery_results.models import TaskResult
from user.models import UserProfile

# others liberys
import json

# Models
from .models import *

# Tasks
from sirna_seeker.tasks import *
from celery.result import AsyncResult



# --------------------------------------------------------------------------------------#

@login_required(redirect_field_name='next',login_url="/accounts/login/")
def dashboard(request):
    user_p = request.user
    user_profile, _ = UserProfile.objects.get_or_create(user=user_p)

    user = request.user
    query_id = request.GET.get('query_id')

    meta_data = {}
    total_sirnas = 0
    sirnas = []

    # Verifica se há um query_id na URL
    if query_id:
        # Tenta buscar a task específica com status 'SUCCESS' e o id do task
        user_task = TaskResult.objects.filter(id=query_id, task_creator=user, status='SUCCESS').first()
        if user_task:
            task_result_string = user_task.result
        else:
            # Se não encontrar o query_id ou se o status não for SUCCESS, redireciona
            return HttpResponseRedirect('sirna/search', {'error': 'Query not found or not successful.'})
    else:
        # Se query_id não for passado, busca a task mais recente com status 'SUCCESS'
        user_tasks = TaskResult.objects.filter(task_creator=user, status='SUCCESS').order_by('-date_created')
        if not user_tasks.exists():
            # Se não houver tasks com sucesso para o usuário, redireciona
            return HttpResponseRedirect('sirna/search', {'error': 'No successful tasks found for this user.'})
        latest_user_task = user_tasks.first()
        task_result_string = latest_user_task.result

    try:
        result_data = json.loads(task_result_string)

        if result_data and isinstance(result_data, list):
            table_data = result_data[0].get('table', [])
            sirnas_verified = result_data[0].get('sirna_verified', [])
            total_sirnas = len(sirnas_verified)

            if len(result_data) > 1 and isinstance(result_data[1], dict):
                meta_data = result_data[1]

            print("Debug - meta_data:", meta_data)

            for index, item in enumerate(table_data):
                siRNA_data = {
                    'index': index + 1,
                    'sequence': item.get("Sequencia'", 'N/A'),
                    'pontuacao': item.get('Pontuacao', 'N/A'),
                    'tm': item.get('TM', 'N/A'),
                    'cg': item.get('CG', 'N/A'),
                    'ag': item.get('AG', 'N/A'),
                    'falhas': item.get('Falhas', []),
                    'posicao': item.get('posicao', 'N/A')
                }
                sirnas.append(siRNA_data)

    except json.JSONDecodeError:
        return render(request, 'resultado_task.html', {'error': 'Error decoding task result data.'})

    paginator = Paginator(sirnas, 10)
    page_number = request.GET.get('page', 1)
    try:
        page_obj = paginator.get_page(page_number)
    except PageNotAnInteger:
        page_obj = paginator.page(1)
    except EmptyPage:
        page_obj = paginator.page(paginator.num_pages)

    MAX_QUERIES = 5
    user_queries = TaskResult.objects.filter(task_creator=user, status='SUCCESS').order_by('-date_created')[:MAX_QUERIES]
    
    # Atualizar títulos das queries
    user_queries_with_titles = []
    for query in user_queries:
        try: 
            query_result_data = json.loads(query.result)
            query_meta_data = {}
            if isinstance(query_result_data, list) and len(query_result_data) > 1 and isinstance(query_result_data[1], dict):
                query_meta_data = query_result_data[1]

            query_title = query_meta_data.get('query_title', '') or f"Unnamed"

        except (json.JSONDecodeError, IndexError, TypeError):
            query_title = f"Unnamed"
        user_queries_with_titles.append({'id': query.id, 'title': query_title, 'date_created': query.date_created})
    
    #output
    context = {
        'sirnas': page_obj.object_list,
        'total_sirnas': paginator.count,
        'page_obj': page_obj,
        'meta_data': meta_data,
        'total_sirnas': total_sirnas,
        'user_queries': user_queries_with_titles,
        'user_profile': {
            'name': user_p.get_full_name(),
            'profission': user_profile.profission,
            'profile_image': user_profile.profile_image,
        }
    }

    return render(request, 'core/core.html', context)
