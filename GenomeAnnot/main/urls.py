from django.urls import include, path
from . import views
from .views import GeneDetailView, GeneValidDetailView, GeneUpdateView, CustomUserLoginView
from django.contrib.auth.views import LoginView, LogoutView
app_name = "main"


urlpatterns = [
    # Home :
    path("", views.home, name="home"),
    
    # Login :
    path('login/', CustomUserLoginView.as_view(), name='login'),
    # path('login/', LoginView.as_view(template_name='main/login.html'), name='login'),

    # path('login/', LoginView.as_view(template_name='login.html'), name='login'),
    # path("logout/", LogoutView.as_view(next_page="login"), name="logout"),
    # path('profile/',  login_required(UserView.as_view()), name='profile'),
    # path('signup/', signup, name='signup'),
    # path("login/", CustomUserLoginView.as_view(), name="custom_user_login"),
    # path("logout/", LogoutView.as_view(next_page="login"), name="logout"),
    
    # Blast and blast with sequence:
    path("blast/", views.blast, name="blast"),
    path("blast/<str:sequence>", views.blast, name="blastseq"),
    
    # Add genome :
    path("addGenome", views.addGenome, name="addGenome"),
    
    # Explore and read gene/genome info :
    path("explore/", views.explore, name="explore"),  # change to listview !!!
    path("explore/genome<int:genome_id>", views.genome, name="genome"),  # change to detailsview !!!
    path("explore/gene<int:gene_id>", GeneDetailView.as_view(), name="gene"),
    
    # Annotate gene info :
    path("annotate/", views.annotate, name="annotate"),  # change to list view !!!
    path("annotate/gene<int:gene_id>",GeneUpdateView.as_view(),name="geneAnnot",),
    
    # Validate gene info :
    path("validate/", views.validate, name="validate"),  # change to list view !!!
    path("validate/gene<int:gene_id>",GeneValidDetailView.as_view(),name="geneValid",),
    
    # Admin
    path("administrator/genome/", views.genomeAdmin, name="genomeAdmin"),
    path("administrator/sequence/", views.sequenceAdmin, name="sequenceAdmin"),
    path("administrator/account/", views.accountAdmin, name="accountAdmin"),
]
