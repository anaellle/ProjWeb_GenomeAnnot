from django.urls import include, path
from . import views
from django.contrib.auth.views import LogoutView
from django.contrib.auth import views as auth_views
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
    ChangePasswordView,
    ResetPasswordView,
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
    
    # Update profile :
    path('profile/',  views.profile, name='profile'),

    # Change password (once logged in) :
    path('password-change/', ChangePasswordView.as_view(), name='password_change'),
    
    # Reset password (if forgotten) :
    path('password-reset/', ResetPasswordView.as_view(), name='password_reset'),
    path('password-reset-confirm/<uidb64>/<token>/',
         auth_views.PasswordResetConfirmView.as_view(template_name='main/password/password_reset_confirm.html'),
         name='password_reset_confirm'),
    path('password-reset-complete/',
         auth_views.PasswordResetCompleteView.as_view(template_name='main/password/password_reset_complete.html'),
         name='password_reset_complete'),
    
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
