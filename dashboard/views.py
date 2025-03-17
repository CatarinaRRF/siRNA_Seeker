# Django libraries
from django.shortcuts import render
from django.http import JsonResponse, HttpResponseRedirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.contrib.auth.decorators import login_required
from django_celery_results.models import TaskResult
from user.models import UserProfile

# other libraries
import json

# Models
from .models import *

# Tasks
from sirna_seeker.tasks import *
from celery.result import AsyncResult
from Bio.Seq import Seq



# --------------------------------------------------------------------------------------#

@login_required(redirect_field_name='next', login_url="/accounts/login/")
def dashboard(request):
    user_p = request.user
    user_profile, _ = UserProfile.objects.get_or_create(user=user_p)
    user = request.user
    query_id = request.GET.get('query_id')

    meta_data = {}
    total_sirnas = 0
    sirnas = []

    if query_id:
        user_task = TaskResult.objects.filter(id=query_id, task_creator=user, status='SUCCESS').first()
        if user_task:
            task_result_string = user_task.result
        else:
            return HttpResponseRedirect('sirna/search', {'error': 'Query not found or not successful.'})
    else:
        user_tasks = TaskResult.objects.filter(task_creator=user, status='SUCCESS').order_by('-date_created')
        if not user_tasks.exists():
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

            for index, item in enumerate(table_data):
                sequence = item.get("Sequencia'", 'N/A')
                
                # Obtendo a sequência complementar
                if sequence != 'N/A':
                    complementary_seq = str(Seq(sequence).complement())  # Obtém a fita complementar
                else:
                    complementary_seq = 'N/A'

                siRNA_data = {
                    'index': index + 1,
                    'sequence': sequence,
                    'complementary_sequence': complementary_seq,
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

def get_all_sequences(request):
    user_p = request.user
    user_tasks = TaskResult.objects.filter(task_creator=user_p, status='SUCCESS').order_by('-date_created')
    
    if not user_tasks.exists():
        return JsonResponse({'error': 'No successful tasks found for this user.'}, status=400)

    latest_user_task = user_tasks.first()
    task_result_string = latest_user_task.result

    try:
        result_data = json.loads(task_result_string)
        all_siRNAs = []

        if result_data and isinstance(result_data, list):
            table_data = result_data[0].get('table', [])
            for index, item in enumerate(table_data):
                sequence = item.get("Sequencia'", 'N/A')
                if sequence != 'N/A':
                    complementary_seq = str(Seq(sequence).complement())
                else:
                    complementary_seq = 'N/A'

                siRNA_data = {
                    'index': index + 1,
                    'sequence': sequence,
                    'complementary_sequence': complementary_seq,
                    'pontuacao': item.get('Pontuacao', 'N/A'),
                    'tm': item.get('TM', 'N/A'),
                    'cg': item.get('CG', 'N/A'),
                    'ag': item.get('AG', 'N/A'),
                    'falhas': item.get('Falhas', []),
                    'posicao': item.get('posicao', 'N/A')
                }
                all_siRNAs.append(siRNA_data)

        return JsonResponse(all_siRNAs, safe=False)

    except json.JSONDecodeError:
        return JsonResponse({'error': 'Error decoding task result data.'}, status=400)
    
@login_required(redirect_field_name='next', login_url="/accounts/login/")
def get_selected_sequences(request):
    user_p = request.user
    selected_indices = request.GET.getlist('indices[]')  # Pega os índices das sequências selecionadas
    
    user_tasks = TaskResult.objects.filter(task_creator=user_p, status='SUCCESS').order_by('-date_created')
    
    if not user_tasks.exists():
        return JsonResponse({'error': 'No successful tasks found for this user.'}, status=400)

    latest_user_task = user_tasks.first()
    task_result_string = latest_user_task.result

    try:
        result_data = json.loads(task_result_string)
        selected_siRNAs = []

        if result_data and isinstance(result_data, list):
            table_data = result_data[0].get('table', [])

            # Processa apenas as sequências que foram selecionadas pelos índices
            for index in selected_indices:
                index = int(index) - 1  # Ajusta para índice baseado em 0
                if 0 <= index < len(table_data):
                    item = table_data[index]
                    sequence = item.get("Sequencia'", 'N/A')

                    if sequence != 'N/A':
                        complementary_seq = str(Seq(sequence).complement())  # Obtém a fita complementar
                    else:
                        complementary_seq = 'N/A'

                    siRNA_data = {
                        'index': index + 1,
                        'sequence': sequence,
                        'complementary_sequence': complementary_seq,
                        'pontuacao': item.get('Pontuacao', 'N/A'),
                        'tm': item.get('TM', 'N/A'),
                        'cg': item.get('CG', 'N/A'),
                        'ag': item.get('AG', 'N/A'),
                        'falhas': item.get('Falhas', []),
                        'posicao': item.get('posicao', 'N/A')
                    }
                    selected_siRNAs.append(siRNA_data)

        return JsonResponse(selected_siRNAs, safe=False)

    except json.JSONDecodeError:
        return JsonResponse({'error': 'Error decoding task result data.'}, status=400)


