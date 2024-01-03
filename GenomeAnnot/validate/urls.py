from django.urls import path
from . import views

app_name = "validate"
urlpatterns = [
    path("", views.index, name="index"),
]
