import GenomeAnnot.wsgi.py
from django.db import models

from django.utils import timezone

from django.contrib.auth.models import AbstractBaseUser, PermissionsMixin
from .managers import CustomUserManager
from django.contrib.auth.models import Permission

from django.utils.translation import gettext_lazy as _

from .managers import CustomUserManager

# Create your models here.


class Genome(models.Model):

    class Status(models.IntegerChoices):
        BLANK = (0, "Not annotated")
        IN_WORK = (1, "In work")
        COMPLETE = (2, "Validated")

    id = models.CharField(max_length=200, primary_key=True)
    species = models.CharField(max_length=200)
    strain = models.CharField(max_length=200, blank=True)
    substrain = models.CharField(max_length=200, blank=True)
    status = models.IntegerField(choices=Status.choices, default=Status.BLANK)

    def __str__(self):
        return self.id


class Chromosome(models.Model):

    class Strand(models.IntegerChoices):
        SENSE = (1, "Sense")
        ANTISENSE = (-1, "Antisense")

    id = models.CharField(max_length=200, primary_key=True)
    chromName = models.CharField(max_length=200)
    startPos = models.IntegerField(default=1)
    endPos = models.IntegerField()
    strand = models.IntegerField(choices=Strand.choices, default=Strand.SENSE)

    idGenome = models.ForeignKey(Genome, on_delete=models.CASCADE)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                "idGenome", "chromName", name="unique_name_of_chromosome_in_genome"
            )
        ]

    def __str__(self):
        return self.id


class ChromosomeSeq(models.Model):

    id = models.CharField(max_length=200, primary_key=True)
    sequence = models.CharField(max_length=int(1e10))

    idChrom = models.ForeignKey(Chromosome, on_delete=models.CASCADE)

    def __str__(self):
        return "SEQ_" + self.id


# class User (models.Model):

#     class Role(models.IntegerChoices):
#         READER = (0, 'Reader')
#         ANNOTATOR = (1, 'Annotator')
#         VALIDATOR = (2, 'Validator')
#         ADMIN = (3, 'Admin')

#     email = models.CharField(max_length=200, primary_key=True)
#     firstName = models.CharField(max_length=50)
#     lastName = models.CharField(max_length=50)
#     password = models.CharField(max_length=15)
#     phoneNumber = models.CharField(max_length=12)
#     role = models.IntegerField(choices=Role.choices, default=Role.READER)
#     lastConnexion = models.DateField()

#     def __str__(self):
#         return self.email


class CustomUser(AbstractBaseUser, PermissionsMixin):

    class Role(models.IntegerChoices):
        READER = (0, "Reader")
        ANNOTATOR = (1, "Annotator")
        VALIDATOR = (2, "Validator")
        ADMIN = (3, "Admin")

    email = models.EmailField(_("email address"), unique=True)

    role = models.IntegerField(choices=Role.choices, default=Role.READER)
    researchCentre = models.CharField(max_length=50, blank=True)
    firstName = models.CharField(max_length=50)
    lastName = models.CharField(max_length=50)
    phoneNumber = models.CharField(max_length=12, blank=True)

    is_staff = models.BooleanField(default=False)
    is_active = models.BooleanField(default=True)
    date_joined = models.DateTimeField(default=timezone.now)

    USERNAME_FIELD = "email"
    REQUIRED_FIELDS = ["firstName", "lastName", "role"]

    objects = CustomUserManager()

    def __str__(self):
        return self.email


class Gene(models.Model):

    class Strand(models.IntegerChoices):
        SENSE = (1, "Sense")
        ANTISENSE = (-1, "Antisense")

    class Status(models.IntegerChoices):
        ASSIGNABLE = (0, "Assignable")
        BEING_ANNOTATED = (1, "Being annotated")
        BEING_CORRECTED = (2, "Being corrected")
        SUBMITTED = (3, "Submitted to a validator")
        VALIDATED = (4, "Validated")

    id = models.CharField(max_length=200, primary_key=True)
    geneName = models.CharField(max_length=200)
    geneSymbol = models.CharField(max_length=100, blank=True)
    geneBiotype = models.CharField(max_length=200, blank=True)
    strand = models.IntegerField(choices=Strand.choices, default=Strand.SENSE)
    startPos = models.IntegerField()
    endPos = models.IntegerField()
    descriptionGene = models.CharField(max_length=1000, blank=True)
    status = models.IntegerField(choices=Status.choices, default=Status.ASSIGNABLE)

    idChrom = models.ForeignKey(Chromosome, on_delete=models.CASCADE)
    emailAnnotator = models.ForeignKey(
        CustomUser,
        related_name="email_annotator",
        on_delete=models.CASCADE,
        blank=True,
        null=True,
    )
    emailValidator = models.ForeignKey(
        CustomUser,
        related_name="email_validator",
        on_delete=models.CASCADE,
        blank=True,
        null=True,
    )

    def __str__(self):
        return self.id


class NucleotidicSeq(models.Model):
    id = models.CharField(max_length=200, primary_key=True)
    sequence = models.CharField(max_length=int(1e8))  # au max : 99 999 999 pb

    idGene = models.ForeignKey(Gene, on_delete=models.CASCADE)

    def __str__(self):
        return "SEQ_" + self.id


class Message(models.Model):

    class TypeOfMessage(models.IntegerChoices):
        AUTO = (0, "Automatic")
        USER = (1, "User")

    # id auto-généré par django
    text = models.CharField(max_length=500)
    date = models.DateTimeField(
        auto_now_add=True
    )  # Automatically set the field to now when the object is first created
    type = models.IntegerField(
        choices=TypeOfMessage.choices, default=TypeOfMessage.AUTO
    )

    emailAuthor = models.ForeignKey(
        CustomUser, on_delete=models.CASCADE, blank=True, null=True
    )
    idGene = models.ForeignKey(Gene, on_delete=models.CASCADE)


class Peptide(models.Model):
    id = models.CharField(max_length=200, primary_key=True)
    transcriptName = models.CharField(max_length=200)
    transcriptBiotype = models.CharField(max_length=200, blank=True)
    descriptionPep = models.CharField(max_length=1000, blank=True)

    idGene = models.ForeignKey(Gene, on_delete=models.CASCADE)

    def __str__(self):
        return self.id


class PeptideSeq(models.Model):
    id = models.CharField(max_length=200, primary_key=True)
    sequence = models.CharField(
        max_length=int(1e6)
    )  # au max : 999 999 acides aminés

    idPeptide = models.ForeignKey(Peptide, on_delete=models.CASCADE)

    def __str__(self):
        return "SEQ_" + self.id


# https://docs.djangoproject.com/en/5.0/ref/models/fields/
# https://docs.djangoproject.com/en/5.0/topics/db/models/
