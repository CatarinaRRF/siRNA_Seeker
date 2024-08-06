from django.shortcuts import render, redirect

from .models import Project, Tag 

def project_list(request):
    tags = Tag.objects.all()

    # Recebe os tipos de projeto selecionados via GET
    selected_types = request.GET.getlist('form-type[]')
    
    # Recebe a universidade selecionada via GET
    university_filter = request.GET.get('university')

    # Filtra os projetos baseado nos filtros selecionados
    projects = Project.objects.all()

    if selected_types:
        projects = projects.filter(tags__id__in=selected_types)

    if university_filter:
        projects = projects.filter(university=university_filter)

    context = {
        'projects': projects,
        'tags': tags,
        'selected_types': selected_types,
        'university_filter': university_filter,
    }

    return render(request, 'pages/projects.html', context)


def error_404(request, exception):
    print('hh')
    return render(request, 'errors/error_404.html', status=404)
 
def error_500(request):
    return render(request, 'errors/error_505.html', status=500)