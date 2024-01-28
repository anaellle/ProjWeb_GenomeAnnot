from django import forms
from .models import Gene, Peptide


class GenePeptideForm(forms.ModelForm):
    class Meta:
        model = Gene
        fields = ["geneName", "geneSymbol", "geneBiotype", "description"]

    peptide_transcriptName = forms.CharField(
        max_length=200, required=False, label="Transcript Name"
    )
    peptide_transcriptBiotype = forms.CharField(
        max_length=200, required=False, label="Transcript Biotype"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields["geneName"].required = False
        self.fields["geneSymbol"].required = False
        self.fields["geneBiotype"].required = False
        self.fields["description"].required = False

    def save(self, commit=True):
        gene = super().save(commit=False)
        peptide_transcriptName = self.cleaned_data.get(
            "peptide_transcriptName"
        )
        peptide_transcriptBiotype = self.cleaned_data.get(
            "peptide_transcriptBiotype"
        )

        if peptide_transcriptName or peptide_transcriptBiotype:
            peptide, created = Peptide.objects.get_or_create(idGene=gene)
            peptide.transcriptName = peptide_transcriptName
            peptide.transcriptBiotype = peptide_transcriptBiotype
            peptide.save()

        if commit:
            gene.save()

        return gene
