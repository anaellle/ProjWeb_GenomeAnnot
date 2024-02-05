import django_filters
from django import forms
from .models import Genome, Gene, CustomUser


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
