from django.shortcuts import render


# Create your views here.
def index(request):
    context = {"active_tab": "blast"}
    return render(request, "blast/main_blast.html", context)
