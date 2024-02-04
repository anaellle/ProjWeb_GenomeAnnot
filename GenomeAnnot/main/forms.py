from django.contrib.auth.forms import UserCreationForm, UserChangeForm
from django import forms
from .models import Gene, Peptide, Message, CustomUser


class CustomUserCreationForm(UserCreationForm):
    """A form for creating new users. Includes all the required
    fields, plus a repeated password.
    """
    email = forms.EmailField(max_length=254)
    
    class Meta:
        model = CustomUser
        fields = ["email",
                  "firstName",
                  "lastName",
                  "researchCentre",
                  "phoneNumber",
                  "role",
                  "password1",
                  "password2",
                  ]


class CustomUserChangeForm(UserChangeForm):
    """A form for updating users. Includes all the fields on
    the user, but replaces the password field with admin's
    disabled password hash display field.
    """
    class Meta:
        model = CustomUser
        fields = '__all__'
        

class GeneUpdateForm(forms.ModelForm):
    class Meta:
        model = Gene
        fields = [
            "geneName",
            "geneSymbol",
            "geneBiotype",
            "descriptionGene",
        ]
        widgets = {
            "descriptionGene": forms.Textarea(attrs={"cols": 100, "rows": 5}),
        }


class PeptideUpdateForm(forms.ModelForm):
    class Meta:
        model = Peptide
        fields = [
            "transcriptName",
            "transcriptBiotype",
            # "descriptionPep",
        ]
        # widgets = {
        #    "descriptionPep": forms.Textarea(attrs={"cols": 100, "rows": 5}),
        # }


class CommentForm(forms.ModelForm):
    class Meta:
        model = Message
        fields = ["text"]
        widgets = {
            "text": forms.Textarea(
                attrs={
                    "cols": 100,
                    "rows": 5,
                    "placeholder": "write comment ...",
                    "required": "true",
                }
            ),
        }

