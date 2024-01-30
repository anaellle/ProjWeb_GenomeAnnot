from django import forms
from .models import Gene, Peptide


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
            "description": forms.Textarea(
                attrs={"cols": 100, "rows": 5}
            ),  # Personnaliser les dimensions de textarea
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
