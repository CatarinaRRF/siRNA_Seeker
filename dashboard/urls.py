from django.urls import path
from .views import dashboard, get_all_sequences, get_selected_sequences

urlpatterns = [
    path('dashboard', dashboard, name='dashboard'),
    path('get_all_sequences/', get_all_sequences, name='get_all_sequences'),
    path('get_selected_sequences/', get_selected_sequences, name='get_selected_sequences'),
]