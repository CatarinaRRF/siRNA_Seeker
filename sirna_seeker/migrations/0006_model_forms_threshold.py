# Generated by Django 5.0.6 on 2025-05-12 18:37

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("sirna_seeker", "0005_delete_extrataskinfo"),
    ]

    operations = [
        migrations.AddField(
            model_name="model_forms",
            name="threshold",
            field=models.IntegerField(blank=True, default=60, null=True),
        ),
    ]
