from django.urls import path
from . import views

app_name = "main"
urlpatterns = [
    path("", views.home, name="home"),
    path("explore/", views.explore, name="explore"),
    path("explore/genome<int:genome_id>", views.genome, name="genome"),
    path("explore/gene<int:gene_id>", views.gene, name="gene"),
    path("annotate/", views.annotate, name="annotate"),
    path("annotate/gene<int:gene_id>", views.geneAnnot, name="geneAnnot"),
    path("validate/", views.validate, name="validate"),
    path("validate/gene<int:gene_id>", views.geneValid, name="geneValid"),
    path("blast/", views.blast, name="blast"),
    path("admin/genome/", views.genomeAdmin, name="genomeAdmin"),
    path("admin/sequence/", views.sequenceAdmin, name="sequenceAdmin"),
    path("admin/account/", views.accountAdmin, name="accountAdmin"),
    path("addGenome", views.addGenome, name="addGenome"),
]
