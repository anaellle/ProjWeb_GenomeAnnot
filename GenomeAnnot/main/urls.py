from django.urls import path
from . import views

app_name = "main"
urlpatterns = [
    path("", views.home, name="home"),
    path("explore/", views.explore, name="explore"),
    path("annotate/", views.annotate, name="annotate"),
    path("validate/", views.validate, name="validate"),
    path("blast/", views.blast, name="blast"),
    path("admin/genome/", views.genomeAdmin, name="genomeAdmin"),
    path("admin/sequence/", views.sequenceAdmin, name="sequenceAdmin"),
    path("admin/account/", views.accountAdmin, name="accountAdmin"),
    path("addGenome", views.addGenome, name="addGenome"),
]
