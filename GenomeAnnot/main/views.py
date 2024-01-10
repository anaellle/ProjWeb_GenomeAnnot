from django.shortcuts import render


# Create your views here.
def explore(request):
    context = {"active_tab": "explore"}
    return render(request, "main/main_explore.html", context)


def annotate(request):
    context = {"active_tab": "annotate"}
    return render(request, "main/main_annotate.html", context)


def validate(request):
    context = {"active_tab": "validate"}
    return render(request, "main/main_validate.html", context)


def blast(request):
    context = {"active_tab": "blast"}
    return render(request, "main/main_blast.html", context)
