from django.shortcuts import render


def home(request):
    context = {"active_tab": "home"}
    return render(request, "main/home.html", context)


def explore(request):
    context = {"active_tab": "explore"}
    return render(request, "main/explore/main_explore.html", context)


def annotate(request):
    context = {"active_tab": "annotate"}
    return render(request, "main/annotate/main_annotate.html", context)


def validate(request):
    context = {"active_tab": "validate"}
    return render(request, "main/validate/main_validate.html", context)


def blast(request):
    context = {"active_tab": "blast"}
    return render(request, "main/blast/main_blast.html", context)


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
    return render(request, "main/gene.html", context)


def geneValid(request, gene_id):
    context = {
        "gene_id": gene_id,
        "genome_id": "56426",
        "active_tab": "validate",
        "role": "validator",
    }  # ex of context (no db for now)
    return render(request, "main/gene.html", context)
