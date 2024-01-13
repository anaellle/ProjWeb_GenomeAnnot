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
