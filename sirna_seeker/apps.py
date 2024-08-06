from django.apps import AppConfig

class SirnaSeekerConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "sirna_seeker"
    
    def ready(self):
        import sirna_seeker.signals
