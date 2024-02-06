from django.contrib import admin

from django.contrib.auth.admin import UserAdmin
from .forms import CustomUserCreationForm, CustomUserUpdateFormAdmin
from .models import CustomUser

# Register your models here.

from .models import *

class CustomUserAdmin(UserAdmin):
    add_form = CustomUserCreationForm
    form = CustomUserUpdateFormAdmin
    model = CustomUser
    list_display = ("email", "role", "is_staff", "is_active",)
    list_filter = ("email", "role", "is_staff", "is_active",)
    fieldsets = (
        ("Profile", {"fields": ("email", "password", "role", "firstName", "lastName", "phoneNumber", "researchCentre", "date_joined", "last_login")}),
        ("Permissions", {"fields": ("is_staff", "is_active", "is_superuser", "groups", "user_permissions")}),
    )
    add_fieldsets = (
        (None, {
            "classes": ("wide",),
            "fields": (
                "email", "password1", "password2", "role",
                "firstName", "lastName", "phoneNumber", "researchCentre",
                "date_joined", "last_login", "is_staff", "is_active",
                "is_superuser", "groups", "user_permissions"
            )}
        ),
    )
    search_fields = ("email",)
    ordering = ("email",)

class MessageAdmin(admin.ModelAdmin):
    readonly_fields = ('date',)
    fields = ("idGene","text","type","emailAuthor","date",)

admin.site.register(Genome)
admin.site.register(Chromosome)
admin.site.register(ChromosomeSeq)
admin.site.register(CustomUser, CustomUserAdmin)
admin.site.register(Gene)
admin.site.register(NucleotidicSeq)
admin.site.register(Message, MessageAdmin)
admin.site.register(Peptide)
admin.site.register(PeptideSeq)
