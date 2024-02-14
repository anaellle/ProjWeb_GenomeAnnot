import django_filters
from django import forms
from django.db.models import Q
from .models import Genome, Gene, Peptide, CustomUser


####################################################################################


class AdminGenomeFilter(django_filters.FilterSet):

    id__contains = django_filters.CharFilter(
        field_name="id",
        lookup_expr="icontains",
        label="Genome ID",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    species__contains = django_filters.CharFilter(
        field_name="species",
        lookup_expr="icontains",
        label="Species",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    strain__contains = django_filters.CharFilter(
        field_name="strain",
        lookup_expr="icontains",
        label="Strain",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    substrain__contains = django_filters.CharFilter(
        field_name="substrain",
        lookup_expr="icontains",
        label="Substrain",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )

    class Meta:
        model = Genome
        fields = {
            "id": ["contains"],
            "species": ["contains"],
            "strain": ["contains"],
            "substrain": ["contains"],
            "status": ["exact"],
        }


####################################################################################


class AdminGeneFilter(django_filters.FilterSet):
    emailAnnotator__isnull = django_filters.BooleanFilter(
        field_name="emailAnnotator",
        lookup_expr="isnull",
        exclude=True,
        label="Annotator assigned ",
    )
    #        widget=forms.TextInput(attrs={"placeholder": "Placeholder Email Annotator"}),

    emailValidator__isnull = django_filters.BooleanFilter(
        field_name="emailValidator",
        lookup_expr="isnull",
        exclude=True,
        label="Validator assigned ",
    )
    idChrom__idGenome__id__icontains = django_filters.CharFilter(
        field_name="idChrom__idGenome__id",
        lookup_expr="icontains",
        label="Genome",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    geneName__contains = django_filters.CharFilter(
        field_name="geneName",
        lookup_expr="icontains",
        label="Gene Name",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    id__contains = django_filters.CharFilter(
        field_name="id",
        lookup_expr="icontains",
        label="Gene ID",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    emailAnnotator__email__icontains = django_filters.CharFilter(
        field_name="emailAnnotator__email",
        lookup_expr="icontains",
        label="Annotator email",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    emailValidator__email__icontains = django_filters.CharFilter(
        field_name="emailValidator__email",
        lookup_expr="icontains",
        label="Validator email",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )

    class Meta:
        model = Gene
        fields = {
            "idChrom__idGenome__id": ["icontains"],
            "id": ["contains"],
            "geneName": ["contains"],
            "status": ["exact"],
            "emailAnnotator__email": ["icontains"],
            "emailValidator__email": ["icontains"],
        }


####################################################################################


class AdminAccountFilter(django_filters.FilterSet):
    firstName__contains = django_filters.CharFilter(
        field_name="firstName",
        lookup_expr="icontains",
        label="First Name",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    lastName__contains = django_filters.CharFilter(
        field_name="lastName",
        lookup_expr="icontains",
        label="Last Name",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    email__contains = django_filters.CharFilter(
        field_name="email",
        lookup_expr="icontains",
        label="Email",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    researchCentre__contains = django_filters.CharFilter(
        field_name="researchCentre",
        lookup_expr="icontains",
        label="ResearchCenter",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    phoneNumber__contains = django_filters.CharFilter(
        field_name="phoneNumber",
        lookup_expr="icontains",
        label="Phone Number",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )

    class Meta:
        model = CustomUser
        fields = {
            "firstName": ["contains"],
            "lastName": ["contains"],
            "email": ["contains"],
            "researchCentre": ["contains"],
            "phoneNumber": ["contains"],
            "role": ["exact"],
        }


class AdminAssignFilter(django_filters.FilterSet):
    firstName__contains = django_filters.CharFilter(
        field_name="firstName",
        lookup_expr="icontains",
        label="First Name",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    lastName__contains = django_filters.CharFilter(
        field_name="lastName",
        lookup_expr="icontains",
        label="Last Name",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    email__contains = django_filters.CharFilter(
        field_name="email",
        lookup_expr="icontains",
        label="Email",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )

    class Meta:
        model = CustomUser
        fields = {
            "firstName": ["contains"],
            "lastName": ["contains"],
            "email": ["contains"],
        }


# ####################################################################################
# ## Annotate
# ####################################################################################


class AnnotateFilter(django_filters.FilterSet):

    # idChrom__idGenome__species = django_filters.ChoiceFilter(
    #     field_name="idChrom__idGenome__species",
    #     label="Species",
    #     choices=[(species, species) for species in Genome.objects.values_list('species', flat=True).distinct()],
    #     widget=forms.Select(attrs={"class": "form-control"})
    # )
    idChrom__idGenome__species__icontains = django_filters.CharFilter(
        field_name="idChrom__idGenome__species",
        lookup_expr="icontains",
        label="Species",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    id__contains = django_filters.CharFilter(
        field_name="id",
        lookup_expr="icontains",
        label="Gene ID",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Gene ID",
                "class": "w-100 rounded border-1",
                "type": "search",
            }
        ),
    )
    geneName__contains = django_filters.CharFilter(
        field_name="geneName",
        lookup_expr="icontains",
        label="Gene Name",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    geneSymbol__contains = django_filters.CharFilter(
        field_name="geneSymbol",
        lookup_expr="icontains",
        label="Gene Symbol",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    sequence_gene = django_filters.CharFilter(
        field_name="nucleotidicseq__sequence",
        lookup_expr="icontains",
        label="Gene Sequence",
        widget=forms.TextInput(
            attrs={
                "placeholder": "",
                "title": "Search for a motif in the gene's nucleotide sequence",
            }
        ),
    )
    sequence_pep = django_filters.CharFilter(
        field_name="peptide__peptideseq__sequence",
        lookup_expr="icontains",
        label="Peptide Sequence",
        widget=forms.TextInput(
            attrs={
                "placeholder": "",
                "title": "Search for a motif in the peptide' sequence",
            }
        ),
    )
    # status__exact = django_filters.ChoiceFilter(
    #     field_name="status",
    #     label="Status",
    #     choices=Gene.Status.choices,
    #     widget=forms.Select(attrs={"class": "form-control"}),
    # )

    notannotated = django_filters.BooleanFilter(
        field_name="status",
        method="filter_notannotated",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    inwork = django_filters.BooleanFilter(
        field_name="status",
        method="filter_inwork",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    review = django_filters.BooleanFilter(
        field_name="status",
        method="filter_review",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    submited = django_filters.BooleanFilter(
        field_name="status",
        method="filter_submited",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    validated = django_filters.BooleanFilter(
        field_name="status",
        method="filter_validated",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )

    def filter_notannotated(self, queryset, name, value):
        if value == False:
            return queryset.exclude(status=Gene.Status.NOT_ANNOTATED)
        return queryset

    def filter_submited(self, queryset, name, value):
        if value == False:
            return queryset.exclude(status=Gene.Status.SUBMITTED)
        return queryset

    def filter_inwork(self, queryset, name, value):
        if value == False:
            return queryset.exclude(status=Gene.Status.BEING_ANNOTATED)
        return queryset

    def filter_review(self, queryset, name, value):
        if value == False:
            return queryset.exclude(status=Gene.Status.BEING_CORRECTED)
        return queryset

    def filter_validated(self, queryset, name, value):
        if value == False:
            return queryset.exclude(status=Gene.Status.VALIDATED)
        return queryset

    # class Meta:
    #     model = Gene
    #     fields = {
    #         "idChrom__idGenome__species":  ["exact"],
    #         "id": ["contains"],
    #         "geneName": ["contains"],
    #         "geneSymbol": ["contains"],
    #     }


# ####################################################################################
# ## Validate
# ####################################################################################


class ValidateFilter(django_filters.FilterSet):

    idChrom__idGenome__species__icontains = django_filters.CharFilter(
        field_name="idChrom__idGenome__species",
        lookup_expr="icontains",
        label="Species",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    id__contains = django_filters.CharFilter(
        field_name="id",
        lookup_expr="icontains",
        label="Gene ID",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Gene ID",
                "class": "w-100 rounded border-1",
                "type": "search",
            }
        ),
    )
    geneName__contains = django_filters.CharFilter(
        field_name="geneName",
        lookup_expr="icontains",
        label="Gene Name",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    geneSymbol__contains = django_filters.CharFilter(
        field_name="geneSymbol",
        lookup_expr="icontains",
        label="Gene Symbol",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    sequence_gene = django_filters.CharFilter(
        field_name="nucleotidicseq__sequence",
        lookup_expr="icontains",
        label="Gene Sequence",
        widget=forms.TextInput(
            attrs={
                "placeholder": "",
                "title": "Search for a motif in the gene's nucleotide sequence",
            }
        ),
    )
    sequence_pep = django_filters.CharFilter(
        field_name="peptide__peptideseq__sequence",
        lookup_expr="icontains",
        label="Peptide Sequence",
        widget=forms.TextInput(
            attrs={
                "placeholder": "",
                "title": "Search for a motif in the peptide' sequence",
            }
        ),
    )

    notannotated = django_filters.BooleanFilter(
        field_name="status",
        method="filter_notannotated",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    tovalidate = django_filters.BooleanFilter(
        field_name="status",
        method="filter_tovalidate",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    validated = django_filters.BooleanFilter(
        field_name="status",
        method="filter_validated",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    # status__exact = django_filters.ChoiceFilter(
    #     field_name="status",
    #     label="Status",
    #     choices=Gene.Status.choices,
    #     widget=forms.Select(attrs={"class": "form-control"}),
    # )

    def filter_notannotated(self, queryset, name, value):
        if value == False:
            excluded_status = (
                Q(status=Gene.Status.BEING_ANNOTATED)
                | Q(status=Gene.Status.BEING_CORRECTED)
                | Q(status=Gene.Status.NOT_ANNOTATED)
            )
            queryset = queryset.exclude(excluded_status)
        return queryset

    def filter_tovalidate(self, queryset, name, value):
        if value == False:
            return queryset.exclude(status=Gene.Status.SUBMITTED)
        return queryset

    def filter_validated(self, queryset, name, value):
        if value == False:
            return queryset.exclude(status=Gene.Status.VALIDATED)
        return queryset


# ####################################################################################
# ## Explore
# ####################################################################################


class ExploreGenomeFilter(django_filters.FilterSet):

    # Filters on Genome
    id__contains = django_filters.CharFilter(
        field_name="id",
        lookup_expr="icontains",
        label="Genome ID",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Genome ID",
                "class": "w-100 rounded border-1",
                "type": "search",
            }
        ),
    )
    species__contains = django_filters.CharFilter(
        field_name="species",
        lookup_expr="icontains",
        label="Species",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    strain__contains = django_filters.CharFilter(
        field_name="strain",
        lookup_expr="icontains",
        label="Strain",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    substrain__contains = django_filters.CharFilter(
        field_name="substrain",
        lookup_expr="icontains",
        label="Substrain",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )

    # Filters on Genome' status
    notannotated = django_filters.BooleanFilter(
        field_name="status",
        method="filter_notannotated",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    inwork = django_filters.BooleanFilter(
        field_name="status",
        method="filter_inwork",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    validated = django_filters.BooleanFilter(
        field_name="status",
        method="filter_validated",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )

    def filter_notannotated(self, queryset, name, value):
        if value == False:
            queryset = queryset.exclude(status=Genome.Status.BLANK)
        return queryset

    def filter_inwork(self, queryset, name, value):
        if value == False:
            queryset = queryset.exclude(status=Genome.Status.IN_WORK)
        return queryset

    def filter_validated(self, queryset, name, value):
        if value == False:
            queryset = queryset.exclude(status=Genome.Status.COMPLETE)
        return queryset

    class Meta:
        model = Genome
        fields = {
            "id": ["contains"],
            "species": ["contains"],
            "strain": ["contains"],
            "substrain": ["contains"],
        }


class ExploreGenePepFilter(django_filters.FilterSet):

    # Filters on Genome
    idChrom__idGenome__id__icontains = django_filters.CharFilter(
        field_name="idChrom__idGenome__id",
        lookup_expr="icontains",
        label="Genome ID",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    # Filters on Chromosome
    idChrom__chromName__icontains = django_filters.CharFilter(
        field_name="idChrom__chromName",
        lookup_expr="icontains",
        label="Chromosome",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    # Filters on Gene
    id__contains = django_filters.CharFilter(
        field_name="id",
        lookup_expr="icontains",
        label="Gene ID",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Gene ID",
                "class": "w-100 rounded border-1",
                "type": "search",
            }
        ),
    )
    geneName__contains = django_filters.CharFilter(
        field_name="geneName",
        lookup_expr="icontains",
        label="Gene Name",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    geneSymbol__contains = django_filters.CharFilter(
        field_name="geneSymbol",
        lookup_expr="icontains",
        label="Gene Symbol",
        widget=forms.TextInput(attrs={"placeholder": ""}),
    )
    motif_sequence_gene = django_filters.CharFilter(
        field_name="nucleotidicseq__sequence",
        lookup_expr="icontains",
        label="DNA motif",
        widget=forms.TextInput(
            attrs={
                "placeholder": "",
                "title": "Search for a motif in the gene's nucleotide sequence",
            }
        ),
    )
    # Filters on Peptide
    id_pep = django_filters.CharFilter(
        field_name="peptide__id",
        lookup_expr="icontains",
        label="Peptide ID",
        widget=forms.TextInput(
            attrs={
                "placeholder": "",
            }
        ),
    )
    name_pep = django_filters.CharFilter(
        field_name="peptide__transcriptName",
        lookup_expr="icontains",
        label="Peptide Name",
        widget=forms.TextInput(
            attrs={
                "placeholder": "",
            }
        ),
    )
    motif_sequence_pep = django_filters.CharFilter(
        field_name="peptide__peptideseq__sequence",
        lookup_expr="icontains",
        label="Peptide Sequence",
        widget=forms.TextInput(
            attrs={
                "placeholder": "",
                "title": "Search for a motif in the peptide' sequence",
            }
        ),
    )
    # Filters on Gene' status
    notannotated = django_filters.BooleanFilter(
        field_name="status",
        method="filter_notannotated",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    inwork = django_filters.BooleanFilter(
        field_name="status",
        method="filter_inwork",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )
    validated = django_filters.BooleanFilter(
        field_name="status",
        method="filter_validated",
        widget=forms.CheckboxInput(
            attrs={"class": "form-check-input", "checked": "checked"}
        ),
    )

    def filter_notannotated(self, queryset, name, value):
        if value == False:
            queryset = queryset.exclude(status=Gene.Status.NOT_ANNOTATED)
        return queryset

    def filter_inwork(self, queryset, name, value):
        if value == False:
            excluded_status = (
                Q(status=Gene.Status.BEING_ANNOTATED)
                | Q(status=Gene.Status.BEING_CORRECTED)
                | Q(status=Gene.Status.SUBMITTED)
            )
            queryset = queryset.exclude(excluded_status)
        return queryset

    def filter_validated(self, queryset, name, value):
        if value == False:
            queryset = queryset.exclude(status=Gene.Status.VALIDATED)
        return queryset
