from django import forms
from django.forms import ModelForm
from .models import model_forms

# Entrada de dados
#----------------------------------------------#

choices = (('reynolds', 'reynolds'), 
           ('ui-tei', "ui-tei"), 
           ('amarzguioui', "amarzguioui")) 

class form_search(forms.ModelForm):

    
    # Configurações siRNA
    #----------------------------------------------#
    sequence = forms.FileField(required=True, widget=forms.FileInput(attrs={'class': 'form-control'})) 
    sequence_tag = forms.CharField(required=False, widget=forms.TextInput(attrs={'class': 'form-control'})) 
    autor = forms.ChoiceField(choices=choices, widget=forms.Select(attrs={'class': 'form-control'})) 
    size = forms.IntegerField(widget=forms.NumberInput(attrs={'class': 'form-control'})) 
    include_tm = forms.ChoiceField(choices=[(True, 'Yes'), (False, 'No')], widget=forms.RadioSelect(attrs={'class': 'form-check-input'}), required=True)
    max_tm = forms.FloatField(initial=21.5, widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.01'}))
    threshold = forms.FloatField(initial=0.6, widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.01'}))

    # Blast
    #----------------------------------------------#
    run_blast = forms.BooleanField(required=False, widget=forms.CheckboxInput(attrs={'class': 'form-check-input'}))
    organism = forms.CharField(initial='human (taxid:9606)', required=False, widget=forms.TextInput(attrs={'class': 'form-control'}))
    database = forms.CharField(initial='refseq_rna', required=False, widget=forms.TextInput(attrs={'class': 'form-control'}))
    identity = forms.FloatField(initial=78, required=False, widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.01'}))
    query_cover = forms.FloatField(initial=78, required=False, widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.01'}))

    class Meta:
        model = model_forms
        fields = ['sequence', 'sequence_tag', 'autor', 
                  'size', 'include_tm', 'max_tm', 'threshold', 'run_blast', 'organism', 
                  'database', 'identity', 'query_cover',]



