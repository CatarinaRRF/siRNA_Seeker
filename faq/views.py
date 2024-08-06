from django.shortcuts import render
from django.contrib.auth.decorators import login_required



# Create your views here.
@login_required(redirect_field_name='next',login_url="/faq")
def teste_view (request):  
    print('you do not change yet')
    return render (request, 'faq.html')