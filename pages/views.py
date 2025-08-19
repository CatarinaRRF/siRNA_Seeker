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
    return render(request, 'errors/error_500.html', status=500)

def authors(request):
    grupos = [
        {
            "id": "grupo1",
            "lider": "Dr. Frederico Pittella Silva",
            "descricao": "Foco em Terapia por RNA de interferência (RNAi) em câncer de mama utilizando nanopartículas. Desenvolvimento de estratégias de entrega eficiente e segura de RNA terapêutico.",
            "membros": [
                "Natália Prado da Silva – Pós-doc",
                "Lívia Mara Silva - Doutoranda",
                "Davi Trombini Aleixo - Doutorando",
                "Estael Cruz Cazarim - Doutoranda",
                "Kezia Cristine Barbosa Ferreira - Doutoranda",
                "Lívia do Nascimento Grossi - Doutoranda",
                "Allana Carvalho Silva - Doutoranda"
            ],
            "foto_lider": "./static/img/frederico.png"
        },
        {
            "id": "grupo2",
            "lider": "Silvia Ligorio Fialho",
            "descricao": "Foco em administração intravítrea de proteínas e dispositivos revestidos por nanofibras para tratamentos oftalmológicos, incluindo degeneração macular relacionada à idade.",
            "membros": [
                "Marina França Dias – Pesquisadora",
                "Lutiana Amaral de Melo – Pesquisadora",
                "Maria Carolina Guerra – Pesquisadora",
                "Deivisson Fagundes – Mestrando",
                "Gabriela Pierazoli – Graduação",
                "Amaury Souza – Graduação"
            ],
            "foto_lider": "https://placehold.co/120x120/cccccc/FFFFFF?text=Foto"
        },
        {
            "id": "grupo3",
            "lider": "Mauro Aparecido de Sousa Xavier",
            "descricao": "Pesquisa em desenvolvimento de vacina recombinante multi-epítopos como controle alternativo do carrapato Rhipicephalus microplus.",
            "membros": [
                "Alessandra Xavier – Colaboradora",
                "Viviane Andrade – Colaboradora",
                "Fábio Ribeiro – Colaborador",
                "Léia Cardoso – Doutoranda",
                "Mayka Rabelo – Colaboradora"
            ],
            "foto_lider": "https://placehold.co/120x120/cccccc/FFFFFF?text=Foto"
        },
        {
            "id": "grupo4",
            "lider": "Luciana Maria Silva Lopes",
            "descricao": "Inativação do gene STAT3 por CRISPR/CAS9 em esferoides 3D enriquecidos de células-tronco do câncer para descoberta de novas terapias oncológicas.",
            "membros": [
                "Milene Pereira Moreira – Pesquisadora",
                "Bianca Nataliene Carvalho de Camargos – Mestranda",
                "Julia Guerra Vieira de Sousa – Mestranda"
            ],
            "foto_lider": "https://placehold.co/120x120/cccccc/FFFFFF?text=Foto"
        },
        {
            "id": "grupo5",
            "lider": "Armando da Silva Cunha Junior",
            "descricao": "Pesquisa em vetores não virais e nanofibras poliméricas para tratamento de distrofias hereditárias de retina, DMRI e glaucoma.",
            "membros": [
                "Carolina Nunes da Silva – Pesquisadora",
                "Ana Carolina Guimaraes Ribeiro – Pesquisadora",
                "Thomas Toshio Inoue – Doutorado",
                "Michelle Gouvêa Gomes – Doutorado",
                "Ana Luiza Faria Silvério – IC",
                "Vitória Luísa Lopes de Souza – IC",
                "Larissa Santos Covre Stefanon – ADC 2A"
            ],
            "foto_lider": "https://placehold.co/120x120/cccccc/FFFFFF?text=Foto"
        },
        {
            "id": "grupo6",
            "lider": "Matheus S. Gomes",
            "descricao": "Projetos em PDAL2 e desenvolvimento de ferramentas web para seleção de siRNAs e miRNAs, integrando bioinformática e biologia molecular.",
            "membros": [
                "Catarina Ribeiro Rezende de Freitas – Desenvolvedora Principal",
                "Valdeir Braz de Paula – Desenvolvedor Principal",
                "Dr. Laurence Rodrigues do Amaral",
                "Dr. Pedro Luiz Lima Bertarini"
            ],
            "foto_lider": "./static/img/matheus.png"
        }
    ]

    return render(request, "god.html", {"grupos": grupos})
