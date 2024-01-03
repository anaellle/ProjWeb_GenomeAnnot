from django.shortcuts import render


# Create your views here.
def index(request):
    context = {"active_tab": "validate"}
    return render(request, "validate/main_validate.html", context)
