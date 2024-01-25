from django.shortcuts import render

# Library required for lauching the Blast API
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import SearchIO


def home(request):
    context = {"active_tab": "home"}
    return render(request, "main/home.html", context)


def explore(request):
    context = {
        "active_tab": "explore",
    }
    if request.method == "GET":
        if "submitsearch" in request.GET:
            # get parameters of search
            searchbar = request.GET.get("searchbar")
            type_res = request.GET.get("res_type")
            context["searchbar"] = searchbar
            context["type_res"] = type_res

            if type_res == "genome":
                strain = request.GET.get("strain")
                species = request.GET.get("species")
                """ status_genome_0 = request.GET.get("0")
                status_genome_1 = request.GET.get("1")
                status_genome_2 = request.GET.get("2") """
                context["strain"] = strain
                context["species"] = species
                """ if "0" in request.GET:
                    context["test"] = "check"
                else:
                    context["test"] = "uncheck"
                context["status_genome_0"] = status_genome_0
                context["status_genome_1"] = status_genome_1
                context["status_genome_2"] = status_genome_2 """

            elif type_res == "gene" or type_res == "prot":
                genome = request.GET.get("genome")
                chrom = request.GET.get("chrom")
                motif = request.GET.get("motif")
                seq = request.GET.get("seq")
                """status_gene_0 = request.GET.get("status0")
                status_gene_123 = request.GET.get("status1")
                status_gene_4 = request.GET.get("status4") """
                context["genome"] = genome
                context["chrom"] = chrom
                context["motif"] = motif
                context["seq"] = seq
                """ context["status_gene_0"] = status_gene_0
                context["status_gene_123"] = status_gene_123
                context["status_gene_4"] = status_gene_4 """

    return render(request, "main/explore/main_explore.html", context)


def annotate(request):
    context = {"active_tab": "annotate"}
    if request.method == "GET":
        if "submitsearch" in request.GET:
            # get parameters of search
            searchbar = request.GET.get("searchbar")
            genome = request.GET.get("genome")
            chrom = request.GET.get("chrom")
            motif_gene = request.GET.get("motif_gene")
            motif_prot = request.GET.get("motif_prot")
            context["searchbar"] = searchbar
            context["genome"] = genome
            context["chrom"] = chrom
            context["motif_gene"] = motif_gene
            context["motif_prot"] = motif_prot

    return render(request, "main/annotate/main_annotate.html", context)


def validate(request):
    context = {"active_tab": "validate"}
    if request.method == "GET":
        if "submitsearch" in request.GET:
            # get parameters of search
            searchbar = request.GET.get("searchbar")
            genome = request.GET.get("genome")
            chrom = request.GET.get("chrom")
            motif_gene = request.GET.get("motif_gene")
            motif_prot = request.GET.get("motif_prot")
            context["searchbar"] = searchbar
            context["genome"] = genome
            context["chrom"] = chrom
            context["motif_gene"] = motif_gene
            context["motif_prot"] = motif_prot
    return render(request, "main/validate/main_validate.html", context)

def blast(request):
    context = {"active_tab": "blast"}
    #return render(request, "main/blast/main_blast.html", context)

    if request.method == 'POST':
        sequence = request.POST['sequence']
        program = request.POST['program']

        db = request.POST['database']
        alignments = request.POST['alignments']
        #Request to ncbi blast api, rajouter gestion des erreurs ensuite
        #try:
        #result_handle = NCBIWWW.qblast(program="blastn", database="nt", sequence=sequence, alignments=5, descriptions=5,format_type="HTML") #Parametres de base pour le moment, rajouter un choix apres
        #blast_results = SearchIO.read(result_handle, "blast-xml") #permet recuperation dans le template pour l'affichage
        #blast_result = result_handle.read()
        #result_handle.close()
        #except Exception as e:
            # Gérer les erreurs, par exemple, en renvoyant un message d'erreur à l'utilisateur
            #return render(request, 'error.html', {'error_message': str(e)})
        
        result_handle = NCBIWWW.qblast(program=program, database=db, sequence=sequence, alignments=alignments, descriptions=50,hitlist_size=5)
        blast_results = SearchIO.read(result_handle, "blast-xml")

        # Traiter les résultats et afficher dans le template
        return render(request, 'main/blast/blast_results.html', {"active_tab": "blast",'results': blast_results})

    return render(request, "main/blast/main_blast.html",context)


def genomeAdmin(request):
    context = {"active_tab": "admin", "active_tab_admin": "genome"}
    return render(request, "main/admin/admin_genome.html", context)


def sequenceAdmin(request):
    context = {"active_tab": "admin", "active_tab_admin": "sequence"}
    return render(request, "main/admin/admin_sequence.html", context)


def accountAdmin(request):
    context = {"active_tab": "admin", "active_tab_admin": "account"}
    return render(request, "main/admin/admin_account.html", context)


def addGenome(request):
    return render(request, "main/addGenome/addGenome.html")


def genome(request, genome_id):  # change to details view later
    context = {
        "genome_id": genome_id,
        "active_tab": "explore",
    }  # ex of context (no db for now)
    return render(request, "main/explore/genome.html", context)


def gene(request, gene_id):  # change to details view later
    context = {
        "gene_id": gene_id,
        "genome_id": "56426",
        "active_tab": "explore",
        "role": "reader",
    }  # ex of context (no db for now)
    return render(request, "main/gene.html", context)


def geneAnnot(request, gene_id):  # change to update view later
    context = {
        "gene_id": gene_id,
        "genome_id": "56426",
        "active_tab": "annotate",
        "role": "annotator",
    }  # ex of context (no db for now)
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


def geneValid(request, gene_id):
    context = {
        "gene_id": gene_id,
        "genome_id": "56426",
        "active_tab": "validate",
        "role": "validator",
    }  # ex of context (no db for now)

    if request.method == "POST":
        # if comment
        if "submit_comment" in request.POST:
            comment = request.POST.get("comment")
        elif "submit_reject" in request.POST:
            ...  # change statut of gene/prot to 2
        elif "submit_validate" in request.POST:
            ...  # change statut to 4

    return render(request, "main/gene.html", context)
