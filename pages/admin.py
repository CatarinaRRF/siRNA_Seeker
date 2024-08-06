from django.contrib import admin
from .models import Project, Tag

# Register your models here.
admin.site.register(Project)
admin.site.register(Tag)