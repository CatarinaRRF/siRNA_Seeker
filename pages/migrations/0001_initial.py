# Generated by Django 5.0.6 on 2024-07-17 15:42

from django.db import migrations, models


class Migration(migrations.Migration):
    initial = True

    dependencies = []

    operations = [
        migrations.CreateModel(
            name="Tag",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("name", models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name="Project",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("title", models.CharField(max_length=200)),
                ("description", models.TextField()),
                ("image_url", models.URLField()),
                ("researcher", models.CharField(max_length=100)),
                ("university", models.CharField(max_length=100)),
                ("location", models.CharField(max_length=100)),
                (
                    "tags",
                    models.ManyToManyField(related_name="projects", to="pages.tag"),
                ),
            ],
        ),
    ]