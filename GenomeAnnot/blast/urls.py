from django.urls import path
from . import views

app_name = "blast"
urlpatterns = [
    path("", views.index, name="index"),
]
