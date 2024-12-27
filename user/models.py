from django.contrib.auth.models import User
from django.db import models

class UserProfile(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE, related_name='profile')
    profission = models.CharField(max_length=255, null=True, blank=True)
    institution = models.CharField(max_length=255, null=True, blank=True)
    profile_image = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return f"Profile of {self.user.username}"
