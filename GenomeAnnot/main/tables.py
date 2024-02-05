# myapp/tables.py
import django_tables2 as tables
from django.urls import reverse
from django.utils.html import format_html
from .models import Genome, Gene, CustomUser, Chromosome


##################################################################
class TableGenome(tables.Table):
    nbgene = tables.Column(
        verbose_name="Number of Genes", empty_values=(), orderable=False
    )
    inwork_percentage = tables.Column(
        verbose_name="In Work Genes (%)", empty_values=(), orderable=False
    )
    valid_percentage = tables.Column(
        verbose_name="Validated Genes (%)", empty_values=(), orderable=False
    )

    class Meta:
        template_name = "django_tables2/table.html"
        model = Genome
        fields = ("id", "species", "strain", "substrain", "status")

    def render_id(self, value):
        url = (
            reverse("main:sequenceAdmin")
            + f"?idChrom__idGenome__id__icontains={value}"
        )
        return format_html(
            '<a href="{}" target="_blank" class="nav-link">{}</a>', url, value
        )

    def render_nbgene(self, record):
        total_genes_count = 0
        chromosomes = Chromosome.objects.filter(idGenome=record)
        for chromosome in chromosomes:
            genes = Gene.objects.filter(idChrom=chromosome)
            total_genes_count += genes.count()
        return total_genes_count

    def render_inwork_percentage(self, record):
        total_genes_count = 0
        inwork_count = 0

        chromosomes = Chromosome.objects.filter(idGenome=record)
        for chromosome in chromosomes:
            genes = Gene.objects.filter(idChrom=chromosome)
            total_genes_count += genes.count()
            inwork_count += genes.exclude(status__in=[0, 4]).count()

        if total_genes_count == 0:
            return None
        else:
            inwork_percentage = (inwork_count / total_genes_count) * 100
            return round(inwork_percentage, 2)

    def render_valid_percentage(self, record):
        total_genes_count = 0
        valid_count = 0

        chromosomes = Chromosome.objects.filter(idGenome=record)
        for chromosome in chromosomes:
            genes = Gene.objects.filter(idChrom=chromosome)
            total_genes_count += genes.count()
            valid_count += genes.filter(status=4).count()

        if total_genes_count == 0:
            return None
        else:
            valid_percentage = (valid_count / total_genes_count) * 100
            return round(valid_percentage, 2)

    def render_status(self, value):
        # include status.html not possible :bug
        if value == "Not annotated":
            string = ' <span class="grey"> \
                <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"\
                  fill="currentColor" class="bi bi-dash" viewBox="0 0 16 16">\
                <path d="M4 8a.5.5 0 0 1 .5-.5h7a.5.5 0 0 1 0 1h-7A.5.5 0 0 1 4 8"/> \
              </svg>\
               Not Annotate  </span>'
        elif value == "In work":
            string = '  <span class="orange"> \
                <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" fill="currentColor" \
                    class="bi bi-exclamation-circle" viewBox="0 0 16 16">\
                <path d="M8 15A7 7 0 1 1 8 1a7 7 0 0 1 0 14m0 1A8 8 0 1 0 8 0a8 8 0 0 0 0 16"/>\
                <path d="M7.002 11a1 1 0 1 1 2 0 1 1 0 0 1-2 0M7.1 4.995a.905.905 0 1 1 1.8 0l-.35 3.507a.552.552 0 0 1-1.1 0z"/>\
              </svg>\
              In Work  </span>'
        elif value == "Validated":
            string = '<span class="green">   \
                <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" fill="currentColor" class="bi bi-check-all" viewBox="0 0 16 16">\
                <path d="M8.97 4.97a.75.75 0 0 1 1.07 1.05l-3.99 4.99a.75.75 0 0 1-1.08.02L2.324 8.384a.75.75 0 1 1 1.06-1.06l2.094 2.093L8.95 4.992zm-.92 5.14.92.92a.75.75 0 0 0 1.079-.02l3.992-4.99a.75.75 0 1 0-1.091-1.028L9.477 9.417l-.485-.486z"/>\
              </svg> \
              Validated   </span>'
        return format_html(string)


##################################################################


class TableGene(tables.Table):

    genome_id = tables.Column(
        verbose_name="Genome ID", accessor="idChrom__idGenome", orderable=True
    )
    choose_user = tables.Column(
        verbose_name="Assign gene", empty_values=(), orderable=False
    )
    emailAnnotator = tables.Column(verbose_name="Annotator")
    emailValidator = tables.Column(verbose_name="Validator")

    class Meta:
        template_name = "django_tables2/table.html"
        model = Gene
        fields = ("id", "geneName", "status", "emailAnnotator", "emailValidator")
        sequence = ("genome_id", "...", "choose_user")

    def render_id(self, record):
        url = reverse("main:gene", kwargs={"gene_id": record.id})
        return format_html(
            '<a href="{}" target="_blank" class="nav-link">{}</a>',
            url,
            record.id,
        )

    def render_choose_user(self, record):
        urlannot = reverse(
            "main:assignAdmin", kwargs={"gene_id": record.id, "role": "1"}
        )
        urlvalid = reverse(
            "main:assignAdmin", kwargs={"gene_id": record.id, "role": "2"}
        )
        return format_html(
            '<a href="{}"  class="nav-link">Annotator</a>\
                <a href="{}"  class="nav-link">Validator</a>',
            urlannot,
            urlvalid,
        )

    def render_genome_id(self, value):

        genome_url = reverse("main:genome", kwargs={"genome_id": value})
        return format_html(
            '<a href="{}" target="_blank" class="nav-link">{}</a>', genome_url, value
        )

    def render_status(self, value):
        # include status.html not possible :bug
        if value == "Assignable":
            string = '<span class="grey">  <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" fill="currentColor" class="bi bi-dash" viewBox="0 0 16 16"><path d="M4 8a.5.5 0 0 1 .5-.5h7a.5.5 0 0 1 0 1h-7A.5.5 0 0 1 4 8"/></svg>Not Annotated  </span>'
        elif value == "Being annotated":
            string = '<span class="orange">   <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" fill="currentColor" class="bi bi-three-dots" viewBox="0 0 16 16"><path d="M3 9.5a1.5 1.5 0 1 1 0-3 1.5 1.5 0 0 1 0 3m5 0a1.5 1.5 0 1 1 0-3 1.5 1.5 0 0 1 0 3m5 0a1.5 1.5 0 1 1 0-3 1.5 1.5 0 0 1 0 3"/></svg>In Work  </span>'
        elif value == "Being corrected":
            string = ' <span class="red">  \
                  <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"\
                      fill="currentColor" class="bi bi-chat-square-dots" viewBox="0 0 16 16"><\
            <path d="M14 1a1 1 0 0 1 1 1v8a1 1 0 0 1-1 1h-2.5a2 2 0 0 0-1.6.8L8 14.333 6.1 11.8a2 2 0 0 0-1.6-.8H2a1 1 0 0 1-1-1V2a1 1 0 0 1 1-1zM2 0a2 2 0 0 0-2 2v8a2 2 0 0 0 2 2h2.5a1 1 0 0 1 .8.4l1.9 2.533a1 1 0 0 0 1.6 0l1.9-2.533a1 1 0 0 1 .8-.4H14a2 2 0 0 0 2-2V2a2 2 0 0 0-2-2z"/><path d="M5 6a1 1 0 1 1-2 0 1 1 0 0 1 2 0m4 0a1 1 0 1 1-2 0 1 1 0 0 1 2 0m4 0a1 1 0 1 1-2 0 1 1 0 0 1 2 0"/></svg>Under Review  </span>'
        elif value == "Submitted to a validator":
            string = '            <span class="green">  \
                <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" fill="currentColor" class="bi bi-check" viewBox="0 0 16 16">\
                    <path d="M10.97 4.97a.75.75 0 0 1 1.07 1.05l-3.99 4.99a.75.75 0 0 1-1.08.02L4.324 8.384a.75.75 0 1 1 1.06-1.06l2.094 2.093 3.473-4.425z"/>\
                </svg>\
                 Submited   </span>'
        elif value == "Validated":
            string = '            <span class="green">  \
                <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" fill="currentColor" class="bi bi-check-all" viewBox="0 0 16 16">\
                    <path d="M8.97 4.97a.75.75 0 0 1 1.07 1.05l-3.99 4.99a.75.75 0 0 1-1.08.02L2.324 8.384a.75.75 0 1 1 1.06-1.06l2.094 2.093L8.95 4.992zm-.92 5.14.92.92a.75.75 0 0 0 1.079-.02l3.992-4.99a.75.75 0 1 0-1.091-1.028L9.477 9.417l-.485-.486z"/>\
                  </svg> \
                 Validated  </span>'
        return format_html(string)

    def render_emailAnnotator(self, value):
        user_url = reverse("main:accountAdmin") + f"?email__contains={value}"
        return format_html('<a href="{}" class="nav-link">{}</a>', user_url, value)

    def render_emailValidator(self, value):
        user_url = reverse("main:accountAdmin") + f"?email__contains={value}"
        return format_html('<a href="{}" class="nav-link">{}</a>', user_url, value)


###################################################################


class TableAccount(tables.Table):

    phoneNumber = tables.Column(verbose_name="Phone", orderable=False)

    class Meta:
        template_name = "django_tables2/table.html"
        model = CustomUser
        fields = (
            "firstName",
            "lastName",
            "email",
            "researchCentre",
            "phoneNumber",
            "role",
            "last_login",
        )


###################################################################


class TableAssignAccount(tables.Table):
    assigned_gene_count = tables.Column(
        verbose_name="Number of Assigned Gene", empty_values=(), orderable=False
    )

    finish_gene_count = tables.Column(
        verbose_name="Work Done (%)", empty_values=(), orderable=False
    )

    assign = tables.Column(
        verbose_name="Assign Gene", empty_values=(), orderable=False
    )

    class Meta:
        template_name = "django_tables2/table.html"
        model = CustomUser
        fields = (
            "firstName",
            "lastName",
            "email",
            "last_login",
        )

    def render_assign(self, record):
        return format_html(
            '<form method="get">\
                        <input class="btn" type="submit" name="{}" value="Choose"></input > \
                            </form>',
            record.email,
        )

    def render_assigned_gene_count(self, record):
        user = CustomUser.objects.get(email=record.email)
        if user.role == 1:
            return Gene.objects.filter(emailAnnotator=user).count()
        elif user.role == 2:
            return Gene.objects.filter(emailValidator=user).count()
        else:
            return None

    def render_finish_gene_count(self, record):
        user = CustomUser.objects.get(email=record.email)
        if user.role == 1:
            gene_count = Gene.objects.filter(emailAnnotator=user).count()
            if gene_count != 0:
                pourcentage = (
                    Gene.objects.filter(
                        emailAnnotator=user, status=Gene.Status.SUBMITTED
                    ).count()
                    / gene_count
                    * 100
                )
                return round(pourcentage, 2)
        elif user.role == 2:
            gene_count = Gene.objects.filter(emailValidator=user).count()
            if gene_count != 0:
                pourcentage = (
                    Gene.objects.filter(
                        emailValidator=user, status=Gene.Status.VALIDATED
                    ).count()
                    / gene_count
                    * 100
                )
                return round(pourcentage, 2)
        else:
            return None
