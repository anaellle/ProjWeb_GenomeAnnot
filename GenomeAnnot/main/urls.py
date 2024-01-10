from django.urls import path
from . import views

app_name = "main"
urlpatterns = [
    path("explore/", views.explore, name="explore"),
    path("annotate/", views.annotate, name="annotate"),
    path("validate/", views.validate, name="validate"),
    path("blast/", views.blast, name="blast"),
]
