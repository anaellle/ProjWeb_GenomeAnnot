from django.urls import path
from . import views
from .views import GeneDetailView, GeneValidDetailView, GeneUpdateView

app_name = "main"


urlpatterns = [
    # Home :
    path("", views.home, name="home"),
    # Blast and blast with sequence:
    path("blast/", views.blast, name="blast"),
    path("blast/<str:sequence>", views.blast, name="blastseq"),
    # Add genome :
    path("addGenome", views.addGenome, name="addGenome"),
    # Explore and read gene/genome info :
    path("explore/", views.explore, name="explore"),  # change to listview !!!
    path(
        "explore/genome<str:genome_id>", views.genome, name="genome"
    ),  # change to detailsview !!!
    path("explore/gene<str:gene_id>", GeneDetailView.as_view(), name="gene"),
    # Annotate gene info :
    path("annotate/", views.annotate, name="annotate"),  # change to list view !!!
    path(
        "annotate/gene<str:gene_id>",
        GeneUpdateView.as_view(),
        name="geneAnnot",
    ),
    # Validate gene info :
    path("validate/", views.validate, name="validate"),  # change to list view !!!
    path(
        "validate/gene<str:gene_id>",
        GeneValidDetailView.as_view(),
        name="geneValid",
    ),
    # Admin
    path("administrator/genome/", views.genomeAdmin, name="genomeAdmin"),
    path("administrator/sequence/", views.sequenceAdmin, name="sequenceAdmin"),
    path("administrator/account/", views.accountAdmin, name="accountAdmin"),
]
