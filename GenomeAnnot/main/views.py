from django.shortcuts import render, redirect, get_object_or_404
from django.views.generic import DetailView, UpdateView, CreateView, View
from django.http import HttpResponseRedirect, Http404, HttpResponse

from django.urls import reverse, resolve
from django.urls.exceptions import Resolver404
from django.utils.html import format_html
from django.core.paginator import Paginator
from django.db.models import Q

from django_tables2 import SingleTableView
from django_tables2.views import SingleTableMixin
from django_tables2.paginators import LazyPaginator
from django_filters.views import FilterView
import plotly.graph_objects as go

from .forms import GeneUpdateForm, PeptideUpdateForm, CommentForm, UploadFileForm
from .tables import TableGenome, TableGene, TableAccount, TableAssignAccount
from .models import Gene, Message, Genome, Chromosome, CustomUser, ChromosomeSeq
from .filters import (
    AdminGenomeFilter,
    AdminGeneFilter,
    AdminAccountFilter,
    AdminAssignFilter,
    ExploreGenomeFilter,
    ExploreGenePepFilter,
    AnnotateFilter,
    ValidateFilter,
)

from .insertion import uploadAndFill

from django.contrib import messages
from django.db.models import Q

# Libraries required for sign up
from django.contrib.auth.mixins import AccessMixin
from django.shortcuts import redirect
from django.views.generic.edit import FormView
from .forms import CustomUserCreationForm
from django.contrib.auth import login

# Libraries required to login
from django.contrib.auth.views import LoginView
from django.contrib import messages
from django.urls import reverse_lazy

# Libraries required for profile
from .forms import CustomUserUpdateForm
from django.contrib.auth.decorators import login_required

# Libraries required to change password (once logged in)
from django.contrib.auth.views import PasswordChangeView

# Libraries required to reset password (if forgotten)
from django.contrib.auth.views import PasswordResetView
from django.contrib.messages.views import SuccessMessageMixin

# Libraries required for lauching the Blast API
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import SearchIO


def get_role(request):
    if request.user.is_authenticated:
        return request.user.role
    else:
        return None


##############################################################################################
######### Home, error 404
##############################################################################################


def home(request):
    context = {"active_tab": "home", "role_user": get_role(request)}
    return render(request, "main/home.html", context)


def custom_404(request, exception):  # only visible if debug set to false
    return render(
        request,
        "main/pageNotFound.html",
        status=404,
        context={"role_user": get_role(request)},
    )


##############################################################################################
######### Sign up
##############################################################################################


class SignUpView(AccessMixin, FormView):
    template_name = "main/signUp.html"
    form_class = CustomUserCreationForm
    success_url = reverse_lazy("main:home")

    def dispatch(self, request, *args, **kwargs):
        if self.request.user.is_authenticated:
            # Redirect authenticated users away from the sign-up page, to the home page
            return redirect(self.get_success_url())
        return super().dispatch(request, *args, **kwargs)

    def get_success_url(self):
        return reverse_lazy("main:home")

    def form_valid(self, form):
        user = form.save()
        if user:
            login(self.request, user)

        return super().form_valid(form)


##############################################################################################
######### Login
##############################################################################################


class CustomUserLoginView(LoginView):
    template_name = "main/login.html"
    redirect_authenticated_user = True

    def get_success_url(self):
        return reverse_lazy("main:home")

    def form_invalid(self, form):
        messages.error(self.request, "Invalid email or password")
        return self.render_to_response(self.get_context_data(form=form))


##############################################################################################
######### Profile
##############################################################################################


@login_required(login_url=reverse_lazy("main:login"))
def profile(request):
    if request.method == "POST":
        user_form = CustomUserUpdateForm(request.POST, instance=request.user)

        if user_form.is_valid():
            user_form.save()
            messages.success(request, "Your profile was successfully updated")
            return redirect(to="main:profile")
    else:
        user_form = CustomUserUpdateForm(instance=request.user)

    return render(
        request,
        "main/profile.html",
        {"role_user": get_role(request), "user_form": user_form},
    )


##############################################################################################
######### Change password (once logged in)
##############################################################################################


class ChangePasswordView(SuccessMessageMixin, PasswordChangeView):
    template_name = "main/password/change_password.html"
    success_message = "Successfully Changed Your Password"
    success_url = reverse_lazy("main:home")


##############################################################################################
######### Reset password (if forgotten) --> non-functional (email not set)
##############################################################################################


class ResetPasswordView(SuccessMessageMixin, PasswordResetView):
    template_name = "main/password/password_reset.html"
    email_template_name = "main/password/password_reset_email.html"
    subject_template_name = "main/password/password_reset_subject.txt"
    success_message = (
        "We've emailed you instructions for setting your password, "
        "if an account exists with the email you entered. You should receive them shortly."
        " If you don't receive an email, "
        "please make sure you've entered the address you registered with, and check your spam folder."
    )
    success_url = reverse_lazy("main:login")


##############################################################################################
######### Search views (explore, annotate and validate)
##############################################################################################


class PaginatedFilterViews(View):
    def get_context_data(self, **kwargs):
        context = super(PaginatedFilterViews, self).get_context_data(**kwargs)
        if self.request.GET:
            querystring = self.request.GET.copy()
            if self.request.GET.get("page"):
                del querystring["page"]
            context["querystring"] = querystring.urlencode()
        return context


# def explore(request):
#     context = {"active_tab": "explore", "role_user": get_role(request)}

#     if request.method == "GET":
#         if "submit_download" in request.GET:
#             ...  # download info gene with gene_id

#         if "submitsearch" in request.GET:
#             # get parameters of search
#             context["searchbar"] = request.GET.get("searchbar")
#             context["type_res"] = request.GET.get("res_type")

#             if context["type_res"] == "genome":
#                 context["strain"] = request.GET.get("strain")
#                 context["species"] = request.GET.get("species")
#                 for status in [
#                     "status0_genome",
#                     "status1_genome",
#                     "status2_genome",
#                 ]:
#                     if status in request.GET:
#                         context[status] = "checked"
#                     else:
#                         context[status] = "unchecked"
#                 # get info about genome
#                 context["genomes_info"] = Genome.objects.all()  # TO DO : filter

#             elif context["type_res"] == "gene" or context["type_res"] == "prot":
#                 # get info filter/search
#                 for info in ["genome", "chrom", "motif", "seq"]:
#                     context[info] = request.GET.get(info)
#                 for status in [
#                     "status0",
#                     "status123",
#                     "status4",
#                 ]:
#                     if status in request.GET:
#                         context[status] = "checked"
#                     else:
#                         context[status] = "unchecked"

#                 # get info about gene/prot
#                 genes = Gene.objects.all()  # TO DO : filter
#                 genes_info = []
#                 for gene in genes:
#                     genome = gene.idChrom.idGenome
#                     peptide = gene.peptide_set.first()
#                     gene_info = {
#                         "gene_id": gene.id,
#                         "gene_name": gene.geneName,
#                         "status": gene.status,
#                         "peptide_id": peptide.id if peptide else None,
#                         "peptide_name": peptide.transcriptName if peptide else None,
#                         "genome_id": genome.id,
#                         "genome_species": genome.species,
#                     }
#                     genes_info.append(gene_info)
#                     context["genes_info"] = genes_info

#     return render(request, "main/explore/main_explore.html", context)


class ExploreGenomeView(PaginatedFilterViews, FilterView):
    model = Genome
    template_name = "main/explore/main_exploreGenome.html"
    paginate_by = 20
    filterset_class = ExploreGenomeFilter

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        role_user = get_role(self.request)
        context["role_user"] = role_user  # pas sûre que ça serve
        context["active_tab"] = "explore"
        context["active_subtab"] = "genome"
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        return queryset.order_by("id")


class ExploreGenePepView(PaginatedFilterViews, FilterView):
    model = Gene
    template_name = "main/explore/main_exploreGenePep.html"
    paginate_by = 20
    filterset_class = ExploreGenePepFilter

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        role_user = get_role(self.request)
        context["role_user"] = role_user  # pas sûre que ça serve
        context["active_tab"] = "explore"
        context["active_subtab"] = "gene"
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        return queryset.order_by("id")


# def annotate(request):
#     context = {"active_tab": "annotate", "role_user": get_role(request)}
#     if context["role_user"] != 1 and context["role_user"] != 3:
#         return redirect("main:home")

#     # get info about gene/prot
#     # filter on annotator
#     genes = Gene.objects.filter(emailAnnotator=request.user)
#     genes_info = []
#     for gene in genes:
#         genome = gene.idChrom.idGenome
#         peptide = gene.peptide_set.first()
#         gene_info = {
#             "gene_id": gene.id,
#             "gene_name": gene.geneName,
#             "status": gene.status,
#             "peptide_id": peptide.id if peptide else None,
#             "peptide_name": peptide.transcriptName if peptide else None,
#             "genome_id": genome.id,
#             "genome_species": genome.species,
#         }
#         genes_info.append(gene_info)
#         context["genes_info"] = genes_info

#     if request.method == "GET":
#         if "submitsearch" in request.GET:
#             # get parameters of search
#             for info in [
#                 "searchbar",
#                 "genome",
#                 "chrom",
#                 "motif_gene",
#                 "motif_prot",
#             ]:
#                 context[info] = request.GET.get(info)
#             for status in [
#                 "status0",
#                 "status1",
#                 "status2",
#                 "status3",
#                 "status4",
#             ]:
#                 if status in request.GET:
#                     context[status] = "checked"
#                 else:
#                     context[status] = "unchecked"

#     return render(request, "main/annotate/main_annotate.html", context)


class AnnotateView(AccessMixin, PaginatedFilterViews, FilterView):
    model = Gene
    template_name = "main/annotate/main_annotate.html"
    paginate_by = 20
    filterset_class = AnnotateFilter

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_authenticated or request.user.role not in (
            CustomUser.Role.ANNOTATOR,
            CustomUser.Role.ADMIN,
        ):
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        role_user = get_role(self.request)
        context["role_user"] = role_user
        context["active_tab"] = "annotate"
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        # Apply the filter to get genes assigned to the current user
        user = self.request.user
        if user:
            assigned_genes = queryset.filter(emailAnnotator=user)
            queryset = queryset.filter(Q(id__in=assigned_genes))
        return queryset.order_by("id")


# def validate(request):
#     context = {"active_tab": "validate", "role_user": get_role(request)}
#     if context["role_user"] != 2 and context["role_user"] != 3:
#         return redirect("main:home")
#     # get info about gene/prot
#     # filter on validator
#     genes = Gene.objects.filter(emailValidator=request.user)
#     genes_info = []
#     for gene in genes:
#         genome = gene.idChrom.idGenome
#         peptide = gene.peptide_set.first()
#         gene_info = {
#             "gene_id": gene.id,
#             "gene_name": gene.geneName,
#             "status": gene.status,
#             "peptide_id": peptide.id if peptide else None,
#             "peptide_name": peptide.transcriptName if peptide else None,
#             "genome_id": genome.id,
#             "genome_species": genome.species,
#         }
#         genes_info.append(gene_info)
#         context["genes_info"] = genes_info

#     if request.method == "GET":
#         if "submitsearch" in request.GET:
#             # get parameters of search
#             for info in [
#                 "searchbar",
#                 "genome",
#                 "chrom",
#                 "motif_gene",
#                 "motif_prot",
#             ]:
#                 context[info] = request.GET.get(info)

#             for status in [
#                 "status012",
#                 "status3",
#                 "status4",
#             ]:
#                 if status in request.GET:
#                     context[status] = "checked"
#                 else:
#                     context[status] = "unchecked"
#     return render(request, "main/validate/main_validate.html", context)


class ValidateView(AccessMixin, PaginatedFilterViews, FilterView):
    model = Gene
    template_name = "main/validate/main_validate.html"
    paginate_by = 20
    filterset_class = ValidateFilter

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_authenticated or request.user.role not in (
            CustomUser.Role.VALIDATOR,
            CustomUser.Role.ADMIN,
        ):
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        role_user = get_role(self.request)
        context["role_user"] = role_user
        context["active_tab"] = "validate"
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        # Apply the filter to get genes assigned to the current user
        user = self.request.user
        if user:
            assigned_genes = queryset.filter(emailValidator=user)
            queryset = queryset.filter(Q(id__in=assigned_genes))
        return queryset.order_by("id")


##############################################################################################
######### Blast
##############################################################################################


import re  # regular expression library


def kind_of_sequence(sequence):
    # Vérifier si la séquence est une séquence d'ADN ou de protéine
    if re.match(r"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]*$", sequence):
        return "nuc"
    elif re.match(
        r"^[ABCDEFGHIKLMNPQRSTUVWXYZabcdefghiklmnpqrstuvwxyz]*$", sequence
    ):
        return "prot"
    else:
        return "pb_seq"


def blast(request, sequence=None):
    context = {
        "active_tab": "blast",
        "role_user": get_role(request),
        "sequence": sequence,
    }

    # return render(request, "main/blast/main_blast.html", context)

    if request.method == "POST":
        sequence = request.POST["sequence"]
        program = request.POST["program"]

        # db = request.POST['database']
        alignments = request.POST["alignments"]
        # Request to ncbi blast api, rajouter gestion des erreurs ensuite

        seq = kind_of_sequence(sequence)
        if seq == "pb_seq":
            context["error_message"] = (
                "Please verify that your query is a protein or a nuc sequence"
            )
            return render(request, "main/blast/error_blast.html", context)
        elif seq == "nuc":
            db = "nt"
            if not (
                program == "blastn" or program == "blastx" or program == "tblastx"
            ):
                context["error_message"] = (
                    "Please choose a programm who works with your type of query (nuc)"
                )
                return render(request, "main/blast/error_blast.html", context)
        elif seq == "prot":
            db = "nr"
            if not (program == "blastp" or program == "tblastn"):
                context["error_message"] = (
                    "Please choose a programm who works with your type of query (prot)"
                )

        try:
            result_handle = NCBIWWW.qblast(
                program=program,
                database=db,
                sequence=sequence,
                alignments=alignments,
                descriptions=50,
                hitlist_size=5,
            )
            blast_results = SearchIO.read(result_handle, "blast-xml")
        except Exception as e:
            # Gérer les erreurs, par exemple, en renvoyant un message d'erreur à l'utilisateur

            context["error_message"] = (
                "No API access, please verify your internet connection"
            )
            return render(request, "main/blast/error_blast.html", context)
        if not blast_results:
            return render(
                request,
                "main/blast/error_blast.html",
                {"error_message": "No results found"},
            )
        context["results"] = blast_results
        return render(request, "main/blast/blast_results.html", context)

    return render(request, "main/blast/main_blast.html", context)


##############################################################################################
######### Add Genome
##############################################################################################


def addGenome(request):
    if not request.user.is_authenticated:
        return redirect("main:home")
    context = {"role_user": get_role(request)}
    if request.method == "POST":
        if "submit_addgenome" in request.POST:
            # get parameters
            form = UploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                genomefile = request.FILES.get("genomefile")
                cdsfile = request.FILES.get("cdsfile")
                peptidefile = request.FILES.get("peptidefile")
                genomeName = genomefile.name.split(".")[0]
                print(genomeName)
                if (genomeName in cdsfile.name) and (genomeName in peptidefile.name):

                    uploadAndFill(genomefile, cdsfile, peptidefile)
                    messages.success(
                        request, "Your files were successfully uploaded"
                    )
                else:
                    messages.error(
                        request, "Your files are not from the same genome"
                    )
                return render(request, "main/addGenome/addGenome.html", context)
            else:
                messages.error(
                    request, "Your files were not uploaded, a problem occured"
                )
    return render(request, "main/addGenome/addGenome.html", context)


##############################################################################################
######### Genome visualization
##############################################################################################


def get_color_status(status):
    if status == 0:
        return "#576b5c"
    elif status in [1, 2, 3]:
        return "#FAB431"
    elif status == 4:
        return "#1CB61C"
    else:
        return "black"


### See genome
class GenomeDetailView(DetailView):  ### TO DO : GENE LIST LINK
    model = Genome
    template_name = "main/explore/genome.html"
    pk_url_kwarg = "genome_id"

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        genome = get_object_or_404(Genome, pk=self.kwargs.get("genome_id"))
        if not request.user.is_authenticated:
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        role_user = get_role(self.request)
        context["genome"] = self.object
        context["active_tab"] = "explore"
        context["active_subtab"] = "genome"
        context["role_user"] = role_user

        ## Figure with plotly :
        fig = go.Figure()
        # get all genes and plot informations
        genes = Gene.objects.filter(idChrom__idGenome=self.object).order_by(
            "startPos"
        )
        # do pagination :
        default_gene_per_page = 30
        gene_per_page = self.request.GET.get("genePerPage", default_gene_per_page)
        gene_per_page = int(gene_per_page)
        paginator = Paginator(genes, gene_per_page)
        page_number = self.request.GET.get("page")
        page_obj = paginator.get_page(page_number)
        context["page_obj"] = page_obj
        context["gene_per_page"] = gene_per_page

        decalage_strand_1 = 0
        decalage_strand_2 = 0
        for gene in page_obj:
            strand = gene.strand
            if strand == 1:
                decalage_strand_1 = 0.3 - decalage_strand_1
                decalage = decalage_strand_1
            elif strand == -1:
                decalage_strand_2 = 0.3 - decalage_strand_2
                decalage = decalage_strand_2
            length = gene.endPos - gene.startPos + 1
            text_pos = f"{gene.startPos}-{gene.endPos}"
            line = go.Scatter(
                x=list(range(gene.startPos, gene.endPos + 1)),
                y=[strand + decalage] * length,
                mode="lines",
                name=gene.id,
                line=dict(color=get_color_status(gene.status), width=6),
                opacity=0.8,
                # hoverinfo="none",
                hovertemplate=text_pos,
            )
            fig.add_trace(line)
            url = reverse("main:gene", kwargs={"gene_id": gene.id})
            link_to_gene = format_html(
                '<a href="{}" target="_blank" style="color: rgb(105, 72, 72); ">{}</a>',
                url,
                gene.id,
            )
            # add link to gene info with annotation text :
            fig.add_annotation(
                x=(gene.startPos + gene.endPos) / 2,
                y=gene.strand + decalage,
                text=link_to_gene,
                ay=-20,
                # showarrow=False,
            )

        # add label
        fig.update_layout(
            yaxis_title="Strand",
            xaxis_title="Position",
            yaxis=dict(range=[-1.3, 1.6], dtick=1),
        )
        # convert to html
        plot_genome = fig.to_html(
            full_html=False, default_height=500, default_width=1500
        )
        context["plot_genome"] = plot_genome

        return context

    def get(self, request, *args, **kwargs):
        if "genelist" in request.GET:
            url = (
                reverse("main:exploreGenePep")
                + f"?idChrom__idGenome__id__icontains={self.get_object()}&notannotated=on&inwork=on&validated=on"
            )
            return HttpResponseRedirect(url)
        return super().get(request, *args, **kwargs)


### Get sequence of genome
class GenomeSeqDetailView(DetailView):
    model = Genome
    template_name = "main/explore/genomeSeq.html"
    pk_url_kwarg = "genome_id"

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        genome = get_object_or_404(Genome, pk=self.kwargs.get("genome_id"))
        if not request.user.is_authenticated:
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        role_user = get_role(self.request)
        context["genome"] = self.object
        context["active_tab"] = "explore"
        context["active_subtab"] = "genome"
        context["role_user"] = role_user
        # get all sequence from chromosomes of genome
        chromosomes = Chromosome.objects.filter(
            idGenome=self.kwargs.get("genome_id")
        )
        sequence = ""
        for chrom in chromosomes:
            chrom_seq = ChromosomeSeq.objects.get(idChrom=chrom)
            sequence += chrom_seq.sequence

        context["sequence"] = sequence
        return context


### Download sequence of Genome
class GenomeSeqDownloadView(View):
    def get(self, request, *args, **kwargs):
        genome_id = kwargs.get("genome_id")
        chromosomes = Chromosome.objects.filter(idGenome=genome_id)
        sequence = ""
        for chrom in chromosomes:
            chrom_seq = ChromosomeSeq.objects.get(idChrom=chrom)
            sequence += (
                ">Chromosome dna:chromosome chromosome:"
                + str(chrom.id)
                + ":Chromosome:"
                + str(chrom.startPos)
                + ":"
                + str(chrom.endPos)
                + ":"
                + str(chrom.strand)
                + " REF\n"
                + str(chrom_seq.sequence)
                + "\n"
            )
        response = HttpResponse(sequence, content_type="text/plain")
        response["Content-Disposition"] = (
            f'attachment; filename="seq_{genome_id}.fa"'
        )
        return response


##############################################################################################
######### Gene reading, annotation and validation
##############################################################################################


#### Function used in view :


# get all related information of a gene
def get_gene_related_info(gene, role):
    chrom = gene.idChrom
    genome = chrom.idGenome
    peptide = gene.peptide_set.first()
    geneseq = gene.nucleotidicseq_set.first()
    peptseq = peptide.peptideseq_set.first() if peptide else None
    if role != 0:
        # if not reader : get the messages associated with gene
        messages = Message.objects.filter(idGene=gene).order_by("date")

    else:
        messages = None
    return {
        "genome": genome,
        "chrom": chrom,
        "gene": gene,
        "geneseq": geneseq,
        "peptide": peptide,
        "peptseq": peptseq,
        "messages": messages,
    }


#### The views of gene


### View to read information of a gene
class GeneDetailView(DetailView):
    model = Gene
    template_name = "main/gene.html"
    pk_url_kwarg = "gene_id"

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        gene = get_object_or_404(Gene, pk=self.kwargs.get("gene_id"))
        if not request.user.is_authenticated:
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # get role of user
        role_user = get_role(self.request)
        # get gene info
        gene = self.object
        gene_info = get_gene_related_info(gene, role_user)
        # merge all info
        context_all = {**context, **gene_info}
        # get accessibility to annotation/validation :
        # we allow admin to annotate and valid in order to test easily our code as admin,
        # but an admin can assign a gene to himself only via admin of django, not interface
        # (so just for us to test)
        if (
            role_user == 1 or role_user == 3
        ) and gene.emailAnnotator == self.request.user:
            context_all["accesAnnot"] = True

        if (
            role_user == 2 or role_user == 3
        ) and gene.emailValidator == self.request.user:
            context_all["accesValid"] = True

        # add info for tabs and role
        context_all["active_tab"] = "explore"
        context_all["active_subtab"] = "gene"
        context_all["role"] = "reader"
        context_all["role_user"] = role_user
        return context_all


class GeneUpdateView(UpdateView):
    model = Gene
    pk_url_kwarg = "gene_id"
    template_name = "main/gene.html"
    form_class = GeneUpdateForm

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        gene = get_object_or_404(Gene, pk=self.kwargs.get("gene_id"))
        if (
            not request.user.is_authenticated
            or request.user.role == 0
            or request.user.role == 2
            or gene.emailAnnotator != self.request.user
        ):
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # url to return after success of a form
    def get_success_url(self):
        return self.request.POST.get(
            "previousAnnot", self.get_context_data()["previous_url"]
        )

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # get role of user
        role_user = get_role(self.request)
        # get gene info
        gene = self.object
        gene_info = get_gene_related_info(gene, role_user)
        # merge all info
        context_all = {**context, **gene_info}
        # add peptide form
        peptide_form = (
            PeptideUpdateForm(instance=context_all["peptide"])
            if context_all["peptide"]
            else None
        )
        context_all["peptide_form"] = peptide_form
        # add info for tabs, role and previous url
        context_all["active_tab"] = "annotate"
        context_all["role"] = "annotator"
        context_all["role_user"] = role_user
        context_all["previous_url"] = self.request.META.get(
            "HTTP_REFERER", reverse("main:annotate")
        )
        return context_all

    # get informations of forms and do actions
    def form_valid(self, form):
        # check user can annotate
        role_user = self.get_context_data()["role_user"]
        gene = self.get_context_data()["gene"]
        if (
            role_user == 1 or role_user == 3
        ) and gene.emailAnnotator == self.request.user:
            # save change of gene and peptide :
            gene = form.save(commit=False)
            peptide = gene.peptide_set.first()
            if peptide:
                peptide_form = PeptideUpdateForm(self.request.POST, instance=peptide)
                if peptide_form.is_valid():
                    peptide_form.save()
            # if save and not annotated, status become 1 (in work) :
            if "submit_save" in self.request.POST:
                # for gene :
                if gene.status == 0:
                    gene.status = 1
                    gene.save()
                # for genome :
                genome = self.object.idChrom.idGenome
                if genome.status == 0:
                    genome.status = 1
                    genome.save()
            # if submited and not just save, status become 3 (submited) :
            if "submit_submit" in self.request.POST:
                user = self.request.user
                gene.status = 3
                gene.save()
                # the status of the genome become "in work" if not already the case
                genome = self.object.idChrom.idGenome
                if genome.status == 0:
                    genome.status = 1
                    genome.save()
                # automatic message of submission :
                message = Message.objects.create(
                    text="Annotation submitted",
                    idGene=gene,
                    emailAuthor=user,
                )
                message.save()

            return super().form_valid(form)
        else:
            return redirect("main:gene", gene_id=gene)


class GeneValidDetailView(DetailView):
    model = Gene
    template_name = "main/gene.html"
    pk_url_kwarg = "gene_id"
    form_class = CommentForm

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        gene = get_object_or_404(Gene, pk=self.kwargs.get("gene_id"))
        if (
            not request.user.is_authenticated
            or request.user.role == 0
            or request.user.role == 1
            or gene.emailValidator != self.request.user
        ):
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # url to return after success of a form
    def get_success_url(self):
        if "submit_comment" not in self.request.POST:
            previous_url = self.request.session.get(
                "previous_url", reverse("main:validate")
            )
            return previous_url
        else:
            return self.request.path

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super(GeneValidDetailView, self).get_context_data(**kwargs)
        role_user = get_role(self.request)
        # get gene info
        gene = self.get_object()
        gene_info = get_gene_related_info(gene, role_user)
        # merge all info
        context_all = {**context, **gene_info}
        # add form comment
        context_all["comment_form"] = self.form_class()
        # add info for tabs and role
        context_all["active_tab"] = "validate"
        context_all["role"] = "validator"
        context_all["role_user"] = role_user
        self.request.session["previous_url"] = self.request.META.get(
            "HTTP_REFERER", reverse("main:validate")
        )
        self.request.session["genome_id"] = context_all["genome"].id
        return context_all

    # get informations of forms and do actions (post form)
    def post(self, request, *args, **kwargs):
        gene = self.get_object()
        if gene.status == 3:  # if gene can be validated/rejected :
            # verification user can validate/reject/comment
            role_user = self.request.user.role
            user = self.request.user
            if (
                role_user == 2 or role_user == 3
            ) and gene.emailValidator == self.request.user:
                # modify gene status and do comment :
                if "submit_comment_reject" in request.POST:
                    form = CommentForm(request.POST)
                    if form.is_valid():
                        # automatic message of rejection :
                        message = Message.objects.create(
                            text="Rejected",
                            idGene=gene,
                            emailAuthor=user,
                        )
                        message.save()
                        # add comment to DB :
                        comment = form.save(commit=False)
                        comment.idGene = gene
                        comment.type = 1
                        comment.emailAuthor = user
                        comment.save()

                        gene.status = 2  # gene in review
                        gene.save()

                elif "submit_validate" in request.POST:
                    gene.status = 4  # gene validated
                    gene.save()
                    # automatic message of validation :
                    message = Message.objects.create(
                        text="Validated",
                        idGene=gene,
                        emailAuthor=user,
                    )
                    message.save()

                    # status of genome is validated if all of gene from genome are validated
                    # if gene added after validation of genome : it is not taken into account
                    # (genome stay validated), but usually all the genome is added at the same time

                    genome_id = self.request.session.get("genome_id", None)
                    genome = Genome.objects.get(id=genome_id)
                    chromosomes = Chromosome.objects.filter(idGenome=genome)
                    all_genes_status_five = all(
                        gene.status == 4
                        for gene in Gene.objects.filter(idChrom__in=chromosomes)
                    )

                    if all_genes_status_five:
                        genome.status = 2
                        genome.save()

                return HttpResponseRedirect(self.get_success_url())
            else:
                return redirect("main:gene", gene_id=gene)


##############################################################################################
######### Administrator
##############################################################################################


class genomeAdmin(SingleTableMixin, FilterView):
    model = Genome
    table_class = TableGenome
    template_name = "main/admin/admin.html"
    paginate_by = 10
    filterset_class = AdminGenomeFilter

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_authenticated or request.user.role != 3:
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["active_tab"] = "admin"
        context["active_tab_admin"] = "genome"
        context["role_user"] = get_role(self.request)
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        return queryset.order_by("id")


class sequenceAdmin(SingleTableMixin, FilterView):
    model = Gene
    table_class = TableGene
    template_name = "main/admin/admin.html"
    paginate_by = 10
    filterset_class = AdminGeneFilter

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_authenticated or request.user.role != 3:
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["active_tab"] = "admin"
        context["active_tab_admin"] = "sequence"
        context["role_user"] = get_role(self.request)

        return context

    # filter by genome
    def get_queryset(self):
        queryset = super().get_queryset()
        genome_id = self.request.GET.get("idChrom__idGenome__id__icontains")
        if genome_id:
            queryset = queryset.filter(idChrom__idGenome__id=genome_id)
        return queryset.order_by("id")

    def post(self, request, *args, **kwargs):

        if "assign" in request.POST:
            # get gene list filtered : one for those without annotator, one for those without validator
            genes_no_annot = AdminGeneFilter(
                request.GET,
                queryset=Gene.objects.filter(Q(emailAnnotator__isnull=True)),
            ).qs
            genes_no_valid = AdminGeneFilter(
                request.GET,
                queryset=Gene.objects.filter(Q(emailValidator__isnull=True)),
            ).qs
            # get all annotators and
            annotators = CustomUser.objects.filter(role=1)
            validators = CustomUser.objects.filter(role=2)
            genes = [genes_no_annot, genes_no_valid]
            users = [annotators, validators]
            for r in range(2):
                i = 0
                user = users[r]
                if len(user) != 0:
                    for gene in genes[r]:
                        email = user[i % len(user)]
                        if r == 0:
                            gene.emailAnnotator = CustomUser.objects.get(email=email)
                        elif r == 1:
                            gene.emailValidator = CustomUser.objects.get(email=email)
                        gene.save()
                        i += 1

        return super().get(request, *args, **kwargs)


class accountAdmin(SingleTableMixin, FilterView):
    model = CustomUser
    table_class = TableAccount
    template_name = "main/admin/admin.html"
    paginate_by = 10
    filterset_class = AdminAccountFilter

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_authenticated or request.user.role != 3:
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["active_tab"] = "admin"
        context["active_tab_admin"] = "account"
        context["role_user"] = get_role(self.request)
        return context

    # filter by email
    def get_queryset(self):
        queryset = super().get_queryset()
        email = self.request.GET.get("email__contains")
        if email:
            queryset = queryset.filter(email=email)
        return queryset.order_by("email")


class accountAssignAdmin(SingleTableMixin, FilterView):
    model = CustomUser
    table_class = TableAssignAccount
    template_name = "main/admin/admin.html"
    paginate_by = 10
    filterset_class = AdminAssignFilter

    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        gene = get_object_or_404(Gene, pk=self.kwargs.get("gene_id"))
        if not request.user.is_authenticated or request.user.role != 3:
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # url to return after success of a form
    def get_success_url(self):
        previous_url = self.request.session.get("previous_url")
        if previous_url:
            return previous_url
        return reverse("main:sequenceAdmin")

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["active_tab"] = "admin"
        context["active_tab_admin"] = "account"
        context["role_user"] = get_role(self.request)
        previous_url = self.request.META.get(
            "HTTP_REFERER", reverse("main:sequenceAdmin")
        )
        try:
            match = resolve(previous_url)
            view_name = match.view_name
            if view_name == "sequenceAdmin":
                # store url only if page sequence
                self.request.session["previous_url"] = previous_url
            else:
                previous_url = reverse("main:sequenceAdmin")

        except Resolver404:
            previous_url = reverse("main:sequenceAdmin")
        return context

    # filter by role
    def get_queryset(self):
        queryset = super().get_queryset()
        role = self.kwargs.get("role")
        if role:
            queryset = queryset.filter(role=role)
        return queryset.order_by("email")

    # get informations of forms and do actions (get form)
    def get(self, request, *args, **kwargs):
        user = None
        for key in request.GET:
            if "@" in key:
                try:
                    user = CustomUser.objects.get(email=key)
                except CustomUser.DoesNotExist:
                    pass
                break

        gene_id = self.kwargs.get("gene_id")
        role = int(self.kwargs.get("role"))
        # assign gene to user
        if gene_id and role and user:
            gene = Gene.objects.get(id=gene_id)
            if role == 1:
                gene.emailAnnotator = user
            elif role == 2:
                gene.emailValidator = user
            gene.save()
            return HttpResponseRedirect(self.get_success_url())

        return super().get(request, *args, **kwargs)
