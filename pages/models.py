from django.db import models

class Project(models.Model):
    title = models.CharField(max_length=200)
    researcher = models.CharField(max_length=100)
    university = models.CharField(max_length=100)
    location = models.CharField(max_length=100)
    tags = models.ManyToManyField('Tag', related_name='projects')

class Tag(models.Model):
    name = models.CharField(max_length=50)

    def __str__(self):
        return self.name