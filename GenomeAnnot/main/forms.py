from django.contrib.auth.forms import UserCreationForm, UserChangeForm
from django import forms
from .models import Gene, Peptide, Message, CustomUser


class CustomUserCreationForm(UserCreationForm):
    """A form for creating new users. Includes all the required
    fields, plus a repeated password.
    """
    email = forms.EmailField(max_length=254,) # mettre à True si on veut empêcher la modification
    firstName = forms.CharField(max_length=50,
                               required=True,
                               widget=forms.TextInput())
    lastName = forms.CharField(max_length=50,
                               required=True,
                               widget=forms.TextInput())
    researchCentre = forms.CharField(max_length=50,
                                     required=False,
                                     widget=forms.TextInput())
    phoneNumber = forms.CharField(max_length=12,
                                  required=False,
                                  widget=forms.TextInput(attrs={'placeholder':'+33...'}))
    role = forms.ChoiceField(required=True,
                             choices=CustomUser.Role.choices,
                             widget=forms.Select(),
                            )

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

class CustomUserUpdateForm(forms.ModelForm):
    """A form for updating users.
    """
    email = forms.EmailField(max_length=254,
                             required=True,
                             widget=forms.TextInput(attrs={'class': 'form-control mb-3', 'placeholder':'ex: jean.dupont@gmail.com...'}),
                             disabled=False,) # mettre à True si on veut empêcher la modification
    firstName = forms.CharField(max_length=50,
                               required=True,
                               widget=forms.TextInput(attrs={'class': 'form-control mb-3', 'placeholder':'ex: Jean...'}))
    lastName = forms.CharField(max_length=50,
                               required=True,
                               widget=forms.TextInput(attrs={'class': 'form-control mb-3', 'placeholder':'ex: Dupont...'}))
    researchCentre = forms.CharField(max_length=50,
                                     required=False,
                                     widget=forms.TextInput(attrs={'class': 'form-control mb-3', 'placeholder':'ex: CNRS...'}))
    phoneNumber = forms.CharField(max_length=12,
                                  required=False,
                                  widget=forms.TextInput(attrs={'class': 'form-control mb-3', 'placeholder':'+33...'}))
    role = forms.ChoiceField(required=True,
                             choices=CustomUser.Role.choices,
                             widget=forms.Select(attrs={'class': 'form-control','title':"Non-editable, contact a staff member if you need to."}),
                             disabled=True,
                            )    
    class Meta:
        model = CustomUser
        fields = ["email",
                  "firstName",
                  "lastName",
                  "researchCentre",
                  "phoneNumber",
                  "role",
                  ]

class CustomUserCreationFormAdmin(UserCreationForm):
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

class CustomUserUpdateFormAdmin(UserChangeForm):
    """A form for updating users designed for django's admin platform. Includes all the fields on
    the user, but replaces the password field with admin's
    disabled password hash display field.
    """
    email = forms.EmailField(max_length=254)
    
    class Meta:
        model = CustomUser
        fields = "__all__"


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
            "descriptionPep",
        ]

        widgets = {
            "descriptionPep": forms.Textarea(attrs={"cols": 100, "rows": 5}),
        }


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
