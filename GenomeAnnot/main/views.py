from django.shortcuts import render, redirect, get_object_or_404
from django.views.generic import DetailView, UpdateView, CreateView
from django.http import HttpResponseRedirect, Http404
from django.urls import reverse, resolve
from django.urls.exceptions import Resolver404

from django_tables2 import SingleTableView
from django_tables2.views import SingleTableMixin
from django_tables2.paginators import LazyPaginator
from django_filters.views import FilterView


from .forms import GeneUpdateForm, PeptideUpdateForm, CommentForm
from .tables import TableGenome, TableGene, TableAccount, TableAssignAccount
from .models import Gene, Message, Genome, Chromosome, CustomUser
from .filters import (
    AdminGenomeFilter,
    AdminGeneFilter,
    AdminAccountFilter,
    AdminAssignFilter,
)

# from fileParser import file_to_dico
# from insertion import addData

# Library required for sign up
from django.views.generic.edit import FormView
from .forms import CustomUserCreationForm
from django.contrib.auth import login 

# Library required for login
from django.contrib.auth.views import LoginView
from django.contrib import messages
from django.urls import reverse_lazy

# Library required for lauching the Blast API
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import SearchIO

role_user = "admin"



##############################################################################################
######### Home, error 404
##############################################################################################


def home(request):
    context = {"active_tab": "home", "role_user": role_user}
    return render(request, "main/home.html", context)


def custom_404(request, exception):  # only visible if debug set to false
    return render(
        request,
        "main/pageNotFound.html",
        status=404,
        context={"role_user": role_user},
    )


##############################################################################################
######### Sign up
##############################################################################################


class SignUpView(FormView):
    template_name = 'main/signUp.html'
    form_class = CustomUserCreationForm
    redirect_authenticated_user = True
    success_url = reverse_lazy('main:home') # ça fonctionne pas, à corriger
    
    def form_valid(self, form):
        user = form.save()
        if user:
            login(self.request, user)
        
        return super(SignUpView, self).form_valid(form)


##############################################################################################
######### Login
##############################################################################################


class CustomUserLoginView(LoginView):
    template_name = "main/login.html"
    redirect_authenticated_user = True
    
    def get_success_url(self):
        return reverse_lazy('main:home') 
    
    def form_invalid(self, form):
        messages.error(self.request,'Invalid email or password')
        return self.render_to_response(self.get_context_data(form=form))


##############################################################################################
######### Search views (explore, annotate and validate)
##############################################################################################


def explore(request):
    context = {"active_tab": "explore", "role_user": role_user}

    if request.method == "GET":
        if "submit_download" in request.GET:
            ...  # download info gene with gene_id

        if "submitsearch" in request.GET:
            # get parameters of search
            context["searchbar"] = request.GET.get("searchbar")
            context["type_res"] = request.GET.get("res_type")

            if context["type_res"] == "genome":
                context["strain"] = request.GET.get("strain")
                context["species"] = request.GET.get("species")
                for status in [
                    "status0_genome",
                    "status1_genome",
                    "status2_genome",
                ]:
                    if status in request.GET:
                        context[status] = "checked"
                    else:
                        context[status] = "unchecked"
                # get info about genome
                context["genomes_info"] = Genome.objects.all()  # TO DO : filter

            elif context["type_res"] == "gene" or context["type_res"] == "prot":
                # get info filter/search
                for info in ["genome", "chrom", "motif", "seq"]:
                    context[info] = request.GET.get(info)
                for status in [
                    "status0",
                    "status123",
                    "status4",
                ]:
                    if status in request.GET:
                        context[status] = "checked"
                    else:
                        context[status] = "unchecked"

                # get info about gene/prot
                genes = Gene.objects.all()  # TO DO : filter
                genes_info = []
                for gene in genes:
                    genome = gene.idChrom.idGenome
                    peptide = gene.peptide_set.first()
                    gene_info = {
                        "gene_id": gene.id,
                        "gene_name": gene.geneName,
                        "status": gene.status,
                        "peptide_id": peptide.id if peptide else None,
                        "peptide_name": peptide.transcriptName if peptide else None,
                        "genome_id": genome.id,
                        "genome_species": genome.species,
                    }
                    genes_info.append(gene_info)
                    context["genes_info"] = genes_info

    return render(request, "main/explore/main_explore.html", context)


def annotate(request):
    context = {"active_tab": "annotate", "role_user": role_user}

    # get info about gene/prot
    genes = (
        Gene.objects.all()
    )  # TO DO : filter (dont celui pr user, avoir seulement seq assignée)
    genes_info = []
    for gene in genes:
        genome = gene.idChrom.idGenome
        peptide = gene.peptide_set.first()
        gene_info = {
            "gene_id": gene.id,
            "gene_name": gene.geneName,
            "status": gene.status,
            "peptide_id": peptide.id if peptide else None,
            "peptide_name": peptide.transcriptName if peptide else None,
            "genome_id": genome.id,
            "genome_species": genome.species,
        }
        genes_info.append(gene_info)
        context["genes_info"] = genes_info

    if request.method == "GET":
        if "submitsearch" in request.GET:
            # get parameters of search
            for info in [
                "searchbar",
                "genome",
                "chrom",
                "motif_gene",
                "motif_prot",
            ]:
                context[info] = request.GET.get(info)
            for status in [
                "status0",
                "status1",
                "status2",
                "status3",
                "status4",
            ]:
                if status in request.GET:
                    context[status] = "checked"
                else:
                    context[status] = "unchecked"

    return render(request, "main/annotate/main_annotate.html", context)


def validate(request):
    context = {"active_tab": "validate", "role_user": role_user}

    # get info about gene/prot
    genes = (
        Gene.objects.all()
    )  # TO DO : filter (dont celui pr user, avoir seulement seq assignée)
    genes_info = []
    for gene in genes:
        genome = gene.idChrom.idGenome
        peptide = gene.peptide_set.first()
        gene_info = {
            "gene_id": gene.id,
            "gene_name": gene.geneName,
            "status": gene.status,
            "peptide_id": peptide.id if peptide else None,
            "peptide_name": peptide.transcriptName if peptide else None,
            "genome_id": genome.id,
            "genome_species": genome.species,
        }
        genes_info.append(gene_info)
        context["genes_info"] = genes_info

    if request.method == "GET":
        if "submitsearch" in request.GET:
            # get parameters of search
            for info in [
                "searchbar",
                "genome",
                "chrom",
                "motif_gene",
                "motif_prot",
            ]:
                context[info] = request.GET.get(info)

            for status in [
                "status012",
                "status3",
                "status4",
            ]:
                if status in request.GET:
                    context[status] = "checked"
                else:
                    context[status] = "unchecked"
    return render(request, "main/validate/main_validate.html", context)


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



def blast(request,sequence=None):
    context = {
        "active_tab": "blast",
        "role_user": role_user,
        "sequence": sequence,
    }
    #return render(request, "main/blast/main_blast.html", context)

    if request.method == 'POST':
        sequence = request.POST['sequence']
        program = request.POST['program']

        #db = request.POST['database']
        alignments = request.POST['alignments']
        #Request to ncbi blast api, rajouter gestion des erreurs ensuite
        
        seq = kind_of_sequence(sequence)
        if (seq == "pb_seq"):
            context["error_message"]="Please verify that your query is a protein or a nuc sequence"
            return render(request, 'main/blast/error_blast.html', context)
        elif(seq == "nuc"):
            db="nt"
            if not(program == "blastn" or program == "blastx" or program =="tblastx"):
                context["error_message"]="Please choose a programm who works with your type of query (nuc)"
                return render(request, 'main/blast/error_blast.html',context)
        elif(seq== "prot"):
            db="nr"
            if not(program == "blastp" or program == "tblastn"):
                context["error_message"]="Please choose a programm who works with your type of query (prot)"
                return render(request, 'main/blast/error_blast.html',context)

        try:
            result_handle = NCBIWWW.qblast(program=program, database=db, sequence=sequence, alignments=alignments, descriptions=50,hitlist_size=5)
            blast_results = SearchIO.read(result_handle, "blast-xml")
        except Exception as e:
            # Gérer les erreurs, par exemple, en renvoyant un message d'erreur à l'utilisateur
            context["error_message"]="No API access, please verify your internet connection"
            return render(request, 'main/blast/error_blast.html', context)
        if not blast_results:
            return render(request, 'main/blast/error_blast.html', {'error_message': 'No results found'})
        context["results"]=blast_results
        return render(request, 'main/blast/blast_results.html',context)

    return render(request, "main/blast/main_blast.html",context)


##############################################################################################
######### Add Genome
##############################################################################################


def addGenome(request):
    context = {"role_user": role_user}
    if request.method == "POST":
        if "submit_addgenome" in request.POST:
            # get parameters
            genomefile = request.POST.get("genomefile")
            cdsfile = request.POST.get("cdsfile")
            peptidefile = request.POST.get("peptidefile")
            # python parser to insert into BD : ...
            # print(genomefile)
            # genomeDict = file_to_dico(genomefile)
            # cdsDict = file_to_dico(cdsfile)
            # pepDict = file_to_dico(peptidefile)
            # if genomeDict==-1 or cdsDict==-1 or pepDict==-1 :
            #     print("File is not a .fa file")
            # else :
            #     addData(genomeDict, cdsDict, pepDict)
    return render(request, "main/addGenome/addGenome.html", context)


##############################################################################################
######### Genome visualization
##############################################################################################


def genome(request, genome_id):  # change to details view later
    context = {
        "genome_id": genome_id,
        "active_tab": "explore",
        "role_user": role_user,
    }  # ex of context (no db for now)
    return render(request, "main/explore/genome.html", context)


##############################################################################################
######### Gene reading, annotation and validation
##############################################################################################


#### Not  view, function used in view :


# get all related information of a gene
def get_gene_related_info(gene):
    chrom = gene.idChrom
    genome = chrom.idGenome
    peptide = gene.peptide_set.first()
    geneseq = gene.nucleotidicseq_set.first()
    peptseq = peptide.peptideseq_set.first() if peptide else None
    if role_user != "reader":
        messages = Message.objects.filter(
            idGene=gene
        )  # TO DO : be sure to order by date !!!
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


# all filed fill ?
def are_all_fields_filled(form):
    for field_name, field in form.fields.items():
        if form[field_name].value() in [None, "", []]:
            return False
    return True


#### 3 views of gene


class GeneDetailView(DetailView):
    model = Gene
    template_name = "main/gene.html"
    pk_url_kwarg = "gene_id"

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # get gene info
        gene = self.object
        gene_info = get_gene_related_info(gene)
        # merge all info
        context_all = {**context, **gene_info}
        # add info for tabs and role
        context_all["active_tab"] = "explore"
        context_all["role"] = "reader"
        context_all["role_user"] = role_user
        return context_all


class GeneUpdateView(UpdateView):
    model = Gene
    pk_url_kwarg = "gene_id"
    template_name = "main/gene.html"
    form_class = GeneUpdateForm

    def get_success_url(self):
        return self.request.POST.get(
            "previousAnnot", self.get_context_data()["previous_url"]
        )

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # get gene info
        gene = self.object
        gene_info = get_gene_related_info(gene)
        # merge all info
        context_all = {**context, **gene_info}
        # add peptide form
        peptide_form = (
            PeptideUpdateForm(instance=context_all["peptide"])
            if context_all["peptide"]
            else None
        )
        context_all["peptide_form"] = peptide_form
        # add info for tabs and role
        context_all["active_tab"] = "annotate"
        context_all["role"] = "annotator"
        context_all["role_user"] = role_user
        context_all["previous_url"] = self.request.META.get(
            "HTTP_REFERER", reverse("main:annotate")
        )
        return context_all

    def form_valid(self, form):
        # TO DO : check user can annotate
        gene = form.save(commit=False)
        peptide = gene.peptide_set.first()
        if peptide:
            peptide_form = PeptideUpdateForm(self.request.POST, instance=peptide)
            if peptide_form.is_valid():
                peptide_form.save()
        # if save and not annotated, status become 1 (in work) :
        if "submit_save" in self.request.POST:
            if gene.status == 0:
                gene.status = 1
                gene.save()
            genome = self.object.idChrom.idGenome
            if genome.status == 0:
                genome.status = 1
                genome.save()
        # if submited and not just save, status become 3 (submited) :
        if "submit_submit" in self.request.POST:
            gene.status = 3
            gene.save()
            # the status of the genome become "in work"
            genome = self.object.idChrom.idGenome
            if genome.status == 0:
                genome.status = 1
                genome.save()
            # automatic message of submission :
            message = Message.objects.create(
                text="Annotation submitted",
                idGene=gene,
                emailAuthor=None,  # to change !!!!
            )
            message.save()

        return super().form_valid(form)


class GeneValidDetailView(DetailView):
    model = Gene
    template_name = "main/gene.html"
    pk_url_kwarg = "gene_id"
    form_class = CommentForm

    def get_success_url(self):
        if "submit_comment" not in self.request.POST:
            previous_url = self.request.session.get(
                "previous_url", reverse("main:validate")
            )
            return previous_url
        else:
            return self.request.path

    def get_context_data(self, **kwargs):
        context = super(GeneValidDetailView, self).get_context_data(**kwargs)
        # get gene info
        gene = self.get_object()
        gene_info = get_gene_related_info(gene)
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

    def post(self, request, *args, **kwargs):
        gene = self.get_object()
        if gene.status == 3:  # if gene can be validated/rejected :
            # TO DO : verification user can validate/reject/comment

            if "submit_comment_reject" in request.POST:
                form = CommentForm(request.POST)
                if form.is_valid():
                    # add comment to DB :
                    comment = form.save(commit=False)
                    comment.idGene = gene
                    comment.type = 1
                    comment.emailAuthor = None  ## change !!
                    comment.save()

                gene.status = 2  # gene in review
                gene.save()
                # automatic message of rejection :
                message = Message.objects.create(
                    text="Rejected",
                    idGene=gene,
                    emailAuthor=None,  # to change !!!!
                )
                message.save()

            elif "submit_validate" in request.POST:
                gene.status = 4  # gene validated
                gene.save()
                # automatic message of validation :
                message = Message.objects.create(
                    text="Validated",
                    idGene=gene,
                    emailAuthor=None,  # to change !!!!
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


##############################################################################################
######### Administrator
##############################################################################################


class genomeAdmin(SingleTableMixin, FilterView):
    model = Genome
    table_class = TableGenome
    template_name = "main/admin/admin.html"
    paginate_by = 10
    paginator_class = LazyPaginator
    filterset_class = AdminGenomeFilter

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["active_tab"] = "admin"
        context["active_tab_admin"] = "genome"
        context["role_user"] = role_user

        return context


class sequenceAdmin(SingleTableMixin, FilterView):
    model = Gene
    table_class = TableGene
    template_name = "main/admin/admin.html"
    paginate_by = 10
    paginator_class = LazyPaginator
    filterset_class = AdminGeneFilter

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["active_tab"] = "admin"
        context["active_tab_admin"] = "sequence"
        context["role_user"] = role_user

        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        genome_id = self.request.GET.get("idChrom__idGenome__id__icontains")
        if genome_id:
            queryset = queryset.filter(idChrom__idGenome__id=genome_id)
        return queryset


class accountAdmin(SingleTableMixin, FilterView):
    model = CustomUser
    table_class = TableAccount
    template_name = "main/admin/admin.html"
    paginate_by = 10
    paginator_class = LazyPaginator
    filterset_class = AdminAccountFilter

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["active_tab"] = "admin"
        context["active_tab_admin"] = "account"
        context["role_user"] = role_user
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        email = self.request.GET.get("email__contains")
        if email:
            queryset = queryset.filter(email=email)
        return queryset


class accountAssignAdmin(SingleTableMixin, FilterView):
    model = CustomUser
    table_class = TableAssignAccount
    template_name = "main/admin/admin.html"
    paginate_by = 10
    paginator_class = LazyPaginator
    filterset_class = AdminAssignFilter

    def get_success_url(self):
        previous_url = self.request.session.get("previous_url")
        if previous_url:
            return previous_url
        return reverse("main:sequenceAdmin")

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["active_tab"] = "admin"
        context["active_tab_admin"] = "account"
        context["role_user"] = role_user
        previous_url = self.request.META.get(
            "HTTP_REFERER", reverse("main:sequenceAdmin")
        )
        try:
            match = resolve(previous_url)
            view_name = match.view_name
            if view_name == "sequenceAdmin":
                # store url only if page sequence
                self.request.session["previous_url"] = previous_url
        except Resolver404:
            previous_url = reverse("main:sequenceAdmin")
        return context

    def get_queryset(self):
        queryset = super().get_queryset()
        role = self.kwargs.get("role")
        if role:
            queryset = queryset.filter(role=role)
        return queryset.order_by("email")

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
        if gene_id and role and user:
            gene = Gene.objects.get(id=gene_id)
            if role == 1:
                gene.emailAnnotator = user
            elif role == 2:
                gene.emailValidator = user
            gene.save()
            return HttpResponseRedirect(self.get_success_url())

        return super().get(request, *args, **kwargs)
