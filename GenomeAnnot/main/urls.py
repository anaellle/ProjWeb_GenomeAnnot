from django.urls import include, path
from . import views
from django.contrib.auth.views import LoginView, LogoutView
from .views import (
    GeneDetailView,
    GeneValidDetailView,
    GeneUpdateView,
    genomeAdmin,
    sequenceAdmin,
    accountAdmin,
    accountAssignAdmin,
    CustomUserLoginView,
    SignUpView,
)

app_name = "main"


urlpatterns = [
    # Home :
    path("", views.home, name="home"),
    
    # Login - Logout :
    path('login/', CustomUserLoginView.as_view(), name='login'),
    path('logout/', LogoutView.as_view(next_page='main:login'), name='logout'),
    
    # Sign up :
    path('signUp/', SignUpView.as_view(), name='signUp'),
    # path('profile/',  login_required(UserView.as_view()), name='profile'),

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
    path("administrator/genome/", genomeAdmin.as_view(), name="genomeAdmin"),
    path("administrator/sequence/", sequenceAdmin.as_view(), name="sequenceAdmin"),
    path("administrator/account/", accountAdmin.as_view(), name="accountAdmin"),
    path(
        "administrator/assign<str:gene_id>/<str:role>",
        accountAssignAdmin.as_view(),
        name="assignAdmin",
    ),
]
