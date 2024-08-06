
# Django lovries
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.contrib.auth.decorators import login_required
from django_celery_results.models import TaskResult

# others liberys
import json

# Models
from .models import *

# Tasks
from sirna_seeker.tasks import *
from celery.result import AsyncResult



# --------------------------------------------------------------------------------------#

@login_required(redirect_field_name='next',login_url="/accounts/login/")
def dashboard (request):  
    user = request.user
    user_tasks = TaskResult.objects.filter(task_creator=user).order_by('-date_created')
    

    if not user_tasks.exists():
        return HttpResponseRedirect('sirna/search', {'error': 'No tasks found for this user.'})

    # Get the latest task for the user
    # <django.db.models.fields.TextField: result>
    latest_user_task = user_tasks.first()
    task_result_string = latest_user_task.result

    try:
        # Assume that the `result` field contains JSON-encoded data
        result_data = json.loads(task_result_string)

        # Extract the siRNA data from the first element in the result list
        if result_data and isinstance(result_data, list):
            table_data = result_data[0].get('table', [])
            sirnas_verified = result_data[0].get('sirna_verified', [])
            total_sirnas = len(sirnas_verified)

            # Extract meta data from the second element (position 1)
            meta_data = result_data[1] if len(result_data) > 1 and isinstance(result_data[1], dict) else {}

        else:
            table_data = []
    
        # Prepare the list of siRNA data for the context
        sirnas = []
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

    # Debugging
    print('Debug result_data:', result_data, "***")
    print('Debug sirnas:', sirnas)
    print('Debug meta_data:', meta_data)

    # Pagination
    paginator = Paginator(sirnas, 10)  # Show 10 results per page
    page_number = request.GET.get('page', 1)  # Get the page number from the query parameters
    try:
        page_obj = paginator.get_page(page_number)
    except PageNotAnInteger:
        # If page is not an integer, deliver first page.
        page_obj = paginator.page(1)
    except EmptyPage:
        # If page is out of range (e.g. 9999), deliver last page of results.
        page_obj = paginator.page(paginator.num_pages)

    context = {
        'sirnas': page_obj.object_list,
        'total_sirnas': paginator.count,
        'page_obj': page_obj,
        'meta_data': meta_data,
        'total_sirnas':total_sirnas,
    }

    return render(request, 'core/core.html', context)

@login_required(redirect_field_name='next',login_url="/accounts/login/")
def teste (request):  
    user = request.user
    user_tasks = TaskResult.objects.filter(task_creator=user).order_by('-date_created')
    

    if not user_tasks.exists():
        return HttpResponseRedirect('sirna/search', {'error': 'No tasks found for this user.'})

    # Get the latest task for the user
    # <django.db.models.fields.TextField: result>
    latest_user_task = user_tasks.first()
    task_result_string = latest_user_task.result

    try:
        # Assume that the `result` field contains JSON-encoded data
        result_data = json.loads(task_result_string)

        # Extract the siRNA data from the first element in the result list
        if result_data and isinstance(result_data, list):
            table_data = result_data[0].get('table', [])
            sirnas_verified = result_data[0].get('sirna_verified', [])
            total_sirnas = len(sirnas_verified)

            # Extract meta data from the second element (position 1)
            meta_data = result_data[1] if len(result_data) > 1 and isinstance(result_data[1], dict) else {}

        else:
            table_data = []
    
        # Prepare the list of siRNA data for the context
        sirnas = []
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

    # Debugging
    print('Debug result_data:', result_data, "***")
    print('Debug sirnas:', sirnas)
    print('Debug meta_data:', meta_data)

    # Pagination
    paginator = Paginator(sirnas, 10)  # Show 10 results per page
    page_number = request.GET.get('page', 1)  # Get the page number from the query parameters
    try:
        page_obj = paginator.get_page(page_number)
    except PageNotAnInteger:
        # If page is not an integer, deliver first page.
        page_obj = paginator.page(1)
    except EmptyPage:
        # If page is out of range (e.g. 9999), deliver last page of results.
        page_obj = paginator.page(paginator.num_pages)

    context = {
        'sirnas': page_obj.object_list,
        'total_sirnas': paginator.count,
        'page_obj': page_obj,
        'meta_data': meta_data,
        'total_sirnas':total_sirnas,
    }

    return render(request, 'core/teste.html', context)