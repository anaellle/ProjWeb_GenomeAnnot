from django.urls import path
from . import views
from .views import GeneDetailView, GeneValidDetailView

app_name = "main"


urlpatterns = [
    path("", views.home, name="home"),
    path("explore/", views.explore, name="explore"),  # change to listview
    path(
        "explore/genome<int:genome_id>", views.genome, name="genome"
    ),  # change to detailsview
    path("explore/gene<int:gene_id>", GeneDetailView.as_view(), name="gene"),
    path("annotate/", views.annotate, name="annotate"),  # change to list view
    path(
        "annotate/gene<int:gene_id>", views.geneAnnot, name="geneAnnot"
    ),  # change to update view
    path("validate/", views.validate, name="validate"),  # change to list view
    path(
        "validate/gene<int:gene_id>",
        GeneValidDetailView.as_view(),
        name="geneValid",
    ),
    path("blast/", views.blast, name="blast"),
    path("blast/<str:sequence>", views.blast, name="blastseq"),
    path("administrator/genome/", views.genomeAdmin, name="genomeAdmin"),
    path("administrator/sequence/", views.sequenceAdmin, name="sequenceAdmin"),
    path("administrator/account/", views.accountAdmin, name="accountAdmin"),
    path("addGenome", views.addGenome, name="addGenome"),
]
