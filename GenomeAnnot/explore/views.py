from django.shortcuts import render

# Create your views here.


def index(request):
    context = {"active_tab": "explore"}
    return render(request, "explore/main_explore.html", context)
