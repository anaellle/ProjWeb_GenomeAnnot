from django.shortcuts import render, redirect, get_object_or_404
from django.views.generic import DetailView
from django.http import HttpResponseRedirect, Http404
from django.urls import reverse

from .models import Gene, Message, Genome


# Library required for lauching the Blast API
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import SearchIO

role_user = "admin"


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
                context[
                    "genomes_info"
                ] = Genome.objects.all()  # TO DO : filter

            elif (
                context["type_res"] == "gene" or context["type_res"] == "prot"
            ):
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
                        "peptide_name": peptide.transcriptName
                        if peptide
                        else None,
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


import re  #regular expression library

def kind_of_sequence(sequence):
    # Vérifier si la séquence est une séquence d'ADN ou de protéine
    if re.match(r'^[ACGTURYKMSWBDHVNacgturykmswbdhvn]*$', sequence):
        return "nuc"
    elif re.match(r'^[ABCDEFGHIKLMNPQRSTUVWXYZabcdefghiklmnpqrstuvwxyz]*$', sequence):
        return "prot"
    else:
        return "pb_seq"


def blast(request):
    context = {"active_tab": "blast"}
    #return render(request, "main/blast/main_blast.html", context)

    if request.method == 'POST':
        sequence = request.POST['sequence']
        program = request.POST['program']

        #db = request.POST['database']
        alignments = request.POST['alignments']
        #Request to ncbi blast api, rajouter gestion des erreurs ensuite
        
        match kind_of_sequence(sequence):
            case "pb_seq":
                return render(request, 'main/blast/error_blast.html', {"active_tab": "blast",'error_message': "Please verify that your query is a protein or a nuc sequence"})
            case "nuc":
                db="nt"
                if not(program == "blastn" or program == "blastx" or program =="tblastx"):
                    return render(request, 'main/blast/error_blast.html', {"active_tab": "blast",'error_message': "Please choose a programm who works with your type of query (nuc)"})
            case "prot":
                db="nr"
                if not(program == "blastp" or program == "tblastn"):
                    return render(request, 'main/blast/error_blast.html', {"active_tab": "blast",'error_message': "Please choose a programm who works with your type of query (prot)"})


        try:
            result_handle = NCBIWWW.qblast(program=program, database=db, sequence=sequence, alignments=alignments, descriptions=50,hitlist_size=5)
            blast_results = SearchIO.read(result_handle, "blast-xml")
        except Exception as e:
            # Gérer les erreurs, par exemple, en renvoyant un message d'erreur à l'utilisateur
            return render(request, 'main/blast/error_blast.html', {"active_tab": "blast",'error_message': "No API access, please verify your internet connection"})
        if not blast_results:
            return render(request, 'main/blast/error_blast.html', {'error_message': 'No results found'})
        
        # Traiter les résultats et afficher dans le template
        return render(request, 'main/blast/blast_results.html', {"active_tab": "blast",'results': blast_results})

    return render(request, "main/blast/main_blast.html",context)



def genomeAdmin(request):
    context = {
        "active_tab": "admin",
        "active_tab_admin": "genome",
        "role_user": role_user,
    }
    return render(request, "main/admin/admin_genome.html", context)


def sequenceAdmin(request):
    context = {
        "active_tab": "admin",
        "active_tab_admin": "sequence",
        "role_user": role_user,
    }
    return render(request, "main/admin/admin_sequence.html", context)


def accountAdmin(request):
    context = {
        "active_tab": "admin",
        "active_tab_admin": "account",
        "role_user": role_user,
    }
    return render(request, "main/admin/admin_account.html", context)


def addGenome(request):
    context = {"role_user": role_user}
    if request.method == "POST":
        if "submit_addgenome" in request.POST:
            # get parameters
            genomefile = request.POST.get("genomefile")
            cdsfile = request.POST.get("cdsfile")
            peptidefile = request.POST.get("peptidefile")
            # python parser to insert into BD : ...
    return render(request, "main/addGenome/addGenome.html", context)


def genome(request, genome_id):  # change to details view later
    context = {
        "genome_id": genome_id,
        "active_tab": "explore",
        "role_user": role_user,
    }  # ex of context (no db for now)
    return render(request, "main/explore/genome.html", context)


class GeneDetailView(DetailView):
    model = Gene
    template_name = "main/gene.html"
    pk_url_kwarg = "gene_id"

    # context to extract from DB
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        gene = self.object
        chrom = gene.idChrom
        genome = chrom.idGenome
        peptide = gene.peptide_set.first()
        geneseq = gene.nucleotidicseq_set.first()
        peptseq = peptide.peptideseq_set.first() if peptide else None

        if role_user == "admin":
            messages = Message.objects.filter(
                idGene=gene
            )  # TO DO : be sure to order by date !!!
            context["messages"] = messages

        context["genome"] = genome
        context["chrom"] = chrom
        context["gene"] = gene
        context["geneseq"] = geneseq
        context["peptide"] = peptide
        context["peptseq"] = peptseq
        context["active_tab"] = "explore"
        context["role"] = "reader"
        context["role_user"] = role_user
        return context


def geneAnnot(request, gene_id):  # change to update view later
    gene = get_object_or_404(Gene, pk=gene_id)
    peptide = gene.peptide_set.first()

    context = {
        "genome": gene.idChrom.idGenome,
        "chrom": gene.idChrom,
        "gene": gene,
        "geneseq": gene.nucleotidicseq_set.first(),
        "peptide": peptide,
        "peptseq": peptide.peptideseq_set.first() if peptide else None,
        "messages": Message.objects.filter(
            idGene=gene
        ),  # TO DO : be sure to order by date !!!,
        "active_tab": "annotate",
        "role": "annotator",
        "role_user": role_user,
    }
    if request.method == "POST":
        if "submit_save" in request.POST or "submit_submit" in request.POST:
            # get annotations
            geneName = request.POST.get("geneName")
            geneSymbol = request.POST.get("geneSymbol")
            geneBiotype = request.POST.get("geneBiotype")
            descriGene = request.POST.get("descriGene")
            transcriptName = request.POST.get("transcriptName")
            transcriptBiotype = request.POST.get("transcriptBiotype")
            descriProt = request.POST.get("descriProt")
            if "submit_submit" in request.POST:
                statut = 3  # if submit : status of gene change to "submit" (3)
            # A FAIRE : sinon pas de changement sauf si premiere soumission 0-> 1 (en fonction valeur du statut)

        # A FAIRE : retour sur page de recherche avec filtre conservé

    return render(request, "main/gene.html", context)


class GeneValidDetailView(DetailView):
    model = Gene
    template_name = "main/gene.html"
    pk_url_kwarg = "gene_id"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        gene = self.object
        peptide = gene.peptide_set.first()

        context["gene"] = gene
        context["chrom"] = gene.idChrom
        context["genome"] = gene.idChrom.idGenome
        context["geneseq"] = gene.nucleotidicseq_set.first()
        context["peptide"] = peptide
        context["peptseq"] = (
            peptide.peptideseq_set.first() if peptide else None
        )
        context["messages"] = Message.objects.filter(
            idGene=gene
        )  # TO DO : be sure to order by date !!!
        context["active_tab"] = "validate"
        context["role"] = "validator"
        context["role_user"] = role_user
        return context

    def post(self, request, *args, **kwargs):
        gene = self.get_object()

        if gene.status == 3:
            # TO DO : verification user can validate/reject/comment

            if "submit_reject" in request.POST:
                gene.status = 2
                gene.save()
                message = Message.objects.create(
                    text="Rejected",
                    idGene=gene,
                    emailAuthor=None,  # to change !!!!
                )
                message.save()
                return HttpResponseRedirect(reverse("validate"))

            elif "submit_validate" in request.POST:
                gene.status = 4
                gene.save()
                message = Message.objects.create(
                    text="Validated",
                    idGene=gene,
                    emailAuthor=None,  # to change !!!!
                )
                message.save()
                return HttpResponseRedirect(reverse("validate"))

            elif "submit_comment" in request.POST:
                comment = request.POST.get("comment")
                message = Message.objects.create(
                    text=comment,
                    idGene=gene,
                    type=1,
                    emailAuthor=None,  # to change !!!!
                )
                message.save()
                return HttpResponseRedirect(reverse("validate"))

        return HttpResponseRedirect(reverse("validate"))
