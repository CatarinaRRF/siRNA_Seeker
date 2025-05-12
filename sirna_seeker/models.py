# Libraries
from django.db import models
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
 
# ----------------------------------------------------------------------------------#
#                                      Forms
# ----------------------------------------------------------------------------------#
choices = (
            ('reynolds', 'reynolds'), 
            ('ui-tei', "ui-tei"), 
            ('amarzguioui', "amarzguioui")
                                            )

class model_forms(models.Model):
    #Keys
    # Usuario
    # Sequencias
    # ------------------------------------------------------------------------------#
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    sequence = models.FileField(max_length=10000, blank=True, null=True)
    sequence_tag = models.CharField(max_length=100, blank=True, null=True)

    # Configurações siRNA
    # ------------------------------------------------------------------------------#
    autor = models.CharField(max_length=50, choices=choices, blank=True, null=True)
    size = models.IntegerField(default=21, blank=True, null=True)
    include_tm = models.BooleanField(default=True)
    max_tm = models.FloatField(blank=True, null=True)
    threshold = models.FloatField(default=0.6, blank=True, null=True)

    # Configurações Blast
    # ------------------------------------------------------------------------------#
    run_blast = models.BooleanField(default=True)
    organism = models.CharField(max_length=50, blank=True, null=True)
    database = models.CharField(max_length=50, blank=True, null=True)
    identity = models.FloatField(blank=True, null=True, default=78)
    query_cover = models.FloatField(blank=True, null=True, default=78)


# ----------------------------------------------------------------------------------#
#                                 Tasks Results
# ------------------------------------------------------------------------------#
TaskResult.add_to_class('task_creator', models.ForeignKey(User, 
                                                          on_delete=models.CASCADE, 
                                                          null=True))