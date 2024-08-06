# Librarys
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.core.cache import cache

# Models
from django_celery_results.models import TaskResult
from django.contrib.auth.models import User

# ------------------------------------------------------- #
#                       Signals                           #
# ------------------------------------------------------- #

@receiver(post_save, sender=TaskResult)
def update_task_creator(sender, instance, **kwargs):
    if instance.task_creator_id is not None:
        return
    
    task_id = instance.task_id
    task_creator_username = cache.get(task_id)
    
    if task_creator_username:
        try:
            task_creator = User.objects.get(username=task_creator_username)
            instance.task_creator = task_creator
            print(f"*** o que foi salvo {task_creator} ***")
            instance.save()
        except User.DoesNotExist:
            pass 