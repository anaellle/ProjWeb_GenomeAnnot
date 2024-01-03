from django.shortcuts import render


# Create your views here.
def index(request):
    context = {"active_tab": "annotate"}
    return render(request, "annotate/main_annotate.html", context)
