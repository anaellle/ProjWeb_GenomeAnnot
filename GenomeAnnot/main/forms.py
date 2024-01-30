from django import forms
from .models import Gene, Peptide, Message


class GeneUpdateForm(forms.ModelForm):
    class Meta:
        model = Gene
        fields = [
            "geneName",
            "geneSymbol",
            "geneBiotype",
            "description",
        ]
        widgets = {
            "description": forms.Textarea(attrs={"cols": 100, "rows": 5}),
        }


class PeptideUpdateForm(forms.ModelForm):
    class Meta:
        model = Peptide
        fields = [
            "transcriptName",
            "transcriptBiotype",
            # "description",
        ]
        # widgets = {
        #    "description": forms.Textarea(attrs={"cols": 100, "rows": 5}),
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
