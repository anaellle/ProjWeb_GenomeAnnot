from django.urls import include, path, reverse_lazy
from . import views
from django.contrib.auth.views import LogoutView
from django.contrib.auth import views as auth_views
from .views import (
    ExploreGenomeView,
    ExploreGenePepView,
    GeneDetailView,
    GeneValidDetailView,
    GeneUpdateView,
    GenomeDetailView,
    GenomeSeqDetailView,
    GenomeSeqDownloadView,
    genomeAdmin,
    sequenceAdmin,
    accountAdmin,
    accountAssignAdmin,
    CustomUserLoginView,
    SignUpView,
    ChangePasswordView,
    ResetPasswordView,
    AnnotateView,
    ValidateView,
)

app_name = "main"


urlpatterns = [
    # Home :
    path("", views.home, name="home"),
    
    # Login - Logout :
    path("login/", CustomUserLoginView.as_view(), name="login"),
    path("logout/", LogoutView.as_view(next_page="main:login"), name="logout"),
    
    # Sign up :
    path("signUp/", SignUpView.as_view(), name="signUp"),
    
    # Update profile :
    path("profile/", views.profile, name="profile"),
    
    # Change password (once logged in) :
    path("password-change/", ChangePasswordView.as_view(), name="password_change"),
    
    # Reset password (if forgotten) :
    path("password-reset/", ResetPasswordView.as_view(), name="password_reset"),
    path(
        "password-reset-confirm/<uidb64>/<token>/",
        auth_views.PasswordResetConfirmView.as_view(
            template_name="main/password/password_reset_confirm.html",
            success_url=reverse_lazy("main:password_reset_complete"),
        ),
        name="password_reset_confirm",
    ),
    path(
        "password-reset-complete/",
        auth_views.PasswordResetCompleteView.as_view(
            template_name="main/password/password_reset_complete.html"
        ),
        name="password_reset_complete",
    ),
    
    # Blast and blast with sequence:
    path("blast/", views.blast, name="blast"),
    path("blast/<str:sequence>", views.blast, name="blastseq"),
    
    # Add genome :
    path("addGenome", views.addGenome, name="addGenome"),
    
    # Explore and read gene/genome info :
    # path("explore/", views.explore, name="explore"),  # change to listview !!!
    path("explore/genome/", ExploreGenomeView.as_view(), name="exploreGenome"),  # CHANGED to list view ##
    path("explore/genepep/", ExploreGenePepView.as_view(), name="exploreGenePep"),  # !! Only accessible via URL (no button on the web page yet)
    path("explore/genome<str:genome_id>", GenomeDetailView.as_view(), name="genome"),
    path(
        "explore/genome<str:genome_id>/sequence",
        GenomeSeqDetailView.as_view(),
        name="genomeSeq",
    ),
    path(
        "explore/downloadGenome<str:genome_id>/sequence",
        GenomeSeqDownloadView.as_view(),
        name="downloadGenomeSeq",
    ),
    path("explore/gene<str:gene_id>", GeneDetailView.as_view(), name="gene"),
    
    # Annotate gene info :
    path("annotate/", AnnotateView.as_view(), name="annotate"),  # CHANGED to list view ##
    path(
        "annotate/gene<str:gene_id>",
        GeneUpdateView.as_view(),
        name="geneAnnot",
    ),
    
    # Validate gene info :
    path("validate/", ValidateView.as_view(), name="validate"),  # CHANGED to list view ##
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
