from django.shortcuts import render, redirect, get_object_or_404
from django.views.generic import DetailView, UpdateView, CreateView, View
from django.http import HttpResponseRedirect, Http404, HttpResponse
from django.urls import reverse, resolve
from django.urls.exceptions import Resolver404
from django.utils.html import format_html
from django.core.paginator import Paginator
from django.db.models import Q
from django.contrib import messages
from urllib.parse import urlparse


# import our codes/models :
from .forms import GeneUpdateForm, PeptideUpdateForm, CommentForm, UploadFileForm
from .tables import TableGenome, TableGene, TableAccount, TableAssignAccount
from .models import Gene, Message, Genome, Chromosome, CustomUser, ChromosomeSeq,NucleotidicSeq,PeptideSeq
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


# Library required for tables
from django_tables2.views import SingleTableMixin

# Library required for filter
from django_filters.views import FilterView

# Libraries required for visualisation
import plotly.graph_objects as go

# Libraries required to sign up
from django.contrib.auth.mixins import AccessMixin
from django.shortcuts import redirect
from django.views.generic.edit import FormView
from .forms import CustomUserCreationForm
from django.contrib.auth import login

# Libraries required to log in
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

# Library required for download in csv
import csv

# Libraries required for lauching the Blast API
import re  # regular expression library
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
######### Reset password (if forgotten) --> reset via terminal (email not set)
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


def get_query_page_changed(self):
    if self.request.GET:
        querystring = self.request.GET.copy()  # get parameters in query
        if self.request.GET.get("page"):  # delete the one about page number
            del querystring["page"]
        return querystring.urlencode()  # encode parameters to url
    return None


class PaginatedFilterViews(View):
    ''' A generic view used to specifically manage the pagination of 
    the results list on the Explore, Annotate and Validate pages '''
    
    def get_context_data(self, **kwargs):
        context = super(PaginatedFilterViews, self).get_context_data(**kwargs)
        context["querystring"] = get_query_page_changed(self)
        return context


class ExploreGenomeView(AccessMixin, PaginatedFilterViews, FilterView):
    model = Genome
    template_name = "main/explore/main_exploreGenome.html"
    paginate_by = 20
    filterset_class = ExploreGenomeFilter

    # Return login page if user not unauthenticated
    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_authenticated :
            return redirect("main:login")
        return super().dispatch(request, *args, **kwargs)
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Pagination for request without submit :
        if not self.request.GET:  # (Wrong way of doing things, but functional :)
            context["querystring"] = (
                "id__contains=&submitsearch=&species__contains=&strain__contains=&substrain__contains=&notannotated=on&inwork=on&validated=on"
            )
        role_user = get_role(self.request)
        context["role_user"] = role_user
        context["active_tab"] = "explore"
        context["active_subtab"] = "genome"
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        return queryset.order_by("id")
    
    def post(self,request,*args,**kwargs):
        ''' Manages the download of the filtered list of genomes '''
        
        if "genome_download" in request.POST:
            genomes = ExploreGenomeFilter(request.GET,queryset=Genome.objects.all()).qs
            sequence = ""
            for genome in genomes:
                chromosomes = Chromosome.objects.filter(idGenome=genome.id)
                for chrom in chromosomes:
                    chrom_seq = ChromosomeSeq.objects.get(idChrom=chrom)
                    # Insert line breaks every 60 characters
                    sequence_with_line_breaks = '\n'.join(str(chrom_seq.sequence)[i:i+60] for i in range(0, len(str(chrom_seq.sequence)), 60))
                    sequence += ('>Genome:'+str(genome.id)+' dna:chromosome chromosome:'+str(chrom.id)+
                                ':Chromosome:'+str(chrom.startPos)+':'+str(chrom.endPos)+':'+str(chrom.strand)+
                                ' REF\n'+
                                sequence_with_line_breaks+'\n')

            response = HttpResponse(sequence, content_type='text/plain')
            response["Content-Disposition"] = 'attachment; filename="genomes.fa"'
            return response

        return redirect("main:exploreGenome")

class ExploreGenePepView(AccessMixin, PaginatedFilterViews, FilterView):
    model = Gene
    template_name = "main/explore/main_exploreGenePep.html"
    paginate_by = 20
    filterset_class = ExploreGenePepFilter

    # Return login page if user not unauthenticated
    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_authenticated :
            return redirect("main:login")
        return super().dispatch(request, *args, **kwargs)
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Pagination for request without submit :
        if not self.request.GET:
            context["querystring"] = (
                "id__contains=&submitsearch=&idChrom__idGenome__id__icontains=&idChrom__chromName__icontains=&notannotated=on&inwork=on&validated=on&geneName__contains=&geneSymbol__contains=&motif_sequence_gene=&id_pep=&name_pep=&motif_sequence_pep="
            )

        role_user = get_role(self.request)
        context["role_user"] = role_user
        context["active_tab"] = "explore"
        context["active_subtab"] = "gene"
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        return queryset.order_by("id")

    def get_success_url(self):
        url = self.request.META.get("HTTP_REFERER")
        print(url)
        if url:
            return HttpResponseRedirect(url)
        return reverse("main:exploreGenePep")

    def post(self,request,*args,**kwargs):
        ''' Manages the download of the filtered list of genes and their associated proteins '''
        
        if "genes_download" in request.POST:
            genes = ExploreGenePepFilter(request.GET,queryset=Gene.objects.all()).qs

            response = HttpResponse(content_type="text/csv")
            response["Content-Disposition"] = 'attachment; filename="genes.csv"'

            writer = csv.writer(response)
            writer.writerow(["Genome","Gene ID", "Gene Name","Gene Symbol","Gene Biotype","Strand","Startpos","EndPos","descriptionGene","Chromosome","GeneSeq","PeptideSeq"])

            for gene in genes:
                writer.writerow([gene.idChrom.idGenome,gene.id, gene.geneName,gene.geneSymbol,gene.geneBiotype,gene.strand,gene.startPos,gene.endPos,gene.descriptionGene,gene.idChrom,NucleotidicSeq.objects.filter(idGene=gene.id).first().sequence,PeptideSeq.objects.filter(idPeptide=gene.id).first().sequence])
            return response

        return redirect("main:exploreGenePep")


class AnnotateView(AccessMixin, PaginatedFilterViews, FilterView):
    model = Gene
    template_name = "main/annotate/main_annotate.html"
    paginate_by = 20
    filterset_class = AnnotateFilter

    # Return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_authenticated :
            return redirect("main:login")
        elif request.user.role not in (
            CustomUser.Role.ANNOTATOR,
            CustomUser.Role.ADMIN,
        ):
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Pagination for request without submit :
        if not self.request.GET:
            context["querystring"] = (
                "id__contains=&submitsearch=&idChrom__idGenome__species__icontains=&geneName__contains=&geneSymbol__contains=&sequence_gene=&sequence_pep=&notannotated=on&inwork=on&review=on&submited=on&validated=on"
            )

        role_user = get_role(self.request)
        context["role_user"] = role_user
        context["active_tab"] = "annotate"
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        
        # Apply the filter to get genes assigned to the current annotator
        user = self.request.user
        if user:
            assigned_genes = queryset.filter(emailAnnotator=user)
            queryset = queryset.filter(Q(id__in=assigned_genes))
        return queryset.order_by("id")


class ValidateView(AccessMixin, PaginatedFilterViews, FilterView):
    model = Gene
    template_name = "main/validate/main_validate.html"
    paginate_by = 20
    filterset_class = ValidateFilter

    # Return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_authenticated :
            return redirect("main:login")
        elif request.user.role not in (
            CustomUser.Role.VALIDATOR,
            CustomUser.Role.ADMIN,
        ):
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Pagination for request without submit :
        if not self.request.GET:
            context["querystring"] = (
                "id__contains=&submitsearch=&idChrom__idGenome__species__icontains=&geneName__contains=&geneSymbol__contains=&sequence_gene=&sequence_pep=&notannotated=on&tovalidate=on&validated=on"
            )

        role_user = get_role(self.request)
        context["role_user"] = role_user
        context["active_tab"] = "validate"
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        
        # Apply the filter to get genes assigned to the current validator
        user = self.request.user
        if user:
            assigned_genes = queryset.filter(emailValidator=user)
            queryset = queryset.filter(Q(id__in=assigned_genes))
        return queryset.order_by("id")


##############################################################################################
######### Blast
##############################################################################################


def kind_of_sequence(sequence):
    ''' Checks whether the sequence is DNA or a protein sequence '''
    
    if re.match(r"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]*$", sequence):
        return "nuc"
    elif re.match(
        r"^[ABCDEFGHIKLMNPQRSTUVWXYZabcdefghiklmnpqrstuvwxyz]*$", sequence
    ):
        return "prot"
    else:
        return "pb_seq"


def blast(request, sequence=None):
    ''' Request to ncbi blast api '''
    
    context = {
        "active_tab": "blast",
        "role_user": get_role(request),
        "sequence": sequence,
    }

    if request.method == "POST":
        sequence = request.POST["sequence"]
        program = request.POST["program"]
        alignments = request.POST["alignments"]
        
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

        # Handles errors, returning an error message to the user
        except Exception as e:
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
    '''
    Allow the user to add information in the database by uploading fasta files on the website.
    The files are put in the parser and then the information are added to the database.
    The 3 files need to have the same name minus _cds and _pep in order to be parsed.
    '''
    if not request.user.is_authenticated:
        return redirect("main:home")
    context = {"role_user": get_role(request)}
    if request.method == "POST":
        if "submit_addgenome" in request.POST:
            # get the form
            form = UploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                # get the files
                genomefile = request.FILES.get("genomefile")
                cdsfile = request.FILES.get("cdsfile")
                peptidefile = request.FILES.get("peptidefile")
                genomeName = genomefile.name.split(".")[0]
                print(genomeName)
                # check if the files names are the same
                if (genomeName in cdsfile.name) and (genomeName in peptidefile.name):
                    # parse the files and add to the database 
                    uploadAndFill(genomefile, cdsfile, peptidefile)
                    messages.success(request, 'Your files were successfully uploaded')
                else :
                    # inform the user that the files names are not corrects
                    messages.error(request, "Your files are not from the same genome")
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
class GenomeDetailView(DetailView):
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

        # create visualisation for each gene
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


### Download sequence of a Genome
class GenomeSeqDownloadView(View):
    
    # return home page if url blocked for this user
    def dispatch(self, request, *args, **kwargs):
        genome = get_object_or_404(Genome, pk=self.kwargs.get("genome_id"))
        if not request.user.is_authenticated:
            return redirect("main:home")
        return super().dispatch(request, *args, **kwargs)

    # Download the sequence of a genome
    def get(self, request, *args, **kwargs):
        genome_id = kwargs.get("genome_id")
        chromosomes = Chromosome.objects.filter(idGenome=genome_id)
        sequence = ""
        for chrom in chromosomes:
            chrom_seq = ChromosomeSeq.objects.get(idChrom=chrom)
            # Insert line breaks every 60 characters
            sequence_with_line_breaks = "\n".join(
                str(chrom_seq.sequence)[i : i + 60]
                for i in range(0, len(str(chrom_seq.sequence)), 60)
            )
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
                + sequence_with_line_breaks
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
        context["query_params"] = get_query_page_changed(self)
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
        context["query_params"] = get_query_page_changed(self)
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

            if not self.request.GET.dict():
                genes_no_annot = Gene.objects.filter(
                    Q(emailAnnotator__isnull=True) & ~Q(status=Gene.Status.VALIDATED)
                )
                genes_no_valid = Gene.objects.filter(
                    Q(emailValidator__isnull=True) & ~Q(status=Gene.Status.VALIDATED)
                )
            else:
                # get gene list filtered : one for those without annotator, one for those without validator
                # exclude gene already validated
                genes_no_annot = AdminGeneFilter(
                    request.GET,
                    queryset=Gene.objects.filter(
                        Q(emailAnnotator__isnull=True)
                        & ~Q(status=Gene.Status.VALIDATED)
                    ),
                ).qs
                genes_no_valid = AdminGeneFilter(
                    request.GET,
                    queryset=Gene.objects.filter(
                        Q(emailValidator__isnull=True)
                        & ~Q(status=Gene.Status.VALIDATED)
                    ),
                ).qs

            # get all annotators and validators
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
                        # assign gene :
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
        context["query_params"] = get_query_page_changed(self)
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
        context["query_params"] = get_query_page_changed(self)
        previous_url = self.request.META.get(
            "HTTP_REFERER", reverse("main:sequenceAdmin")
        )
        try:
            # not a good way, but check if view is adminSequence doesn't work, so check url string insteed :
            if "administrator/sequence" in previous_url:
                # store url only if page sequence
                print(previous_url)
                self.request.session["previous_url"] = previous_url
            else:
                self.request.session["previous_url"] = reverse("main:sequenceAdmin")
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
