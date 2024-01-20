from django.db import models

from django.db import models
from django.utils import timezone
from django.urls import reverse
from django.contrib import admin

# Create your models here.

class Genome (models.Model):
    
    class Status(models.IntegerChoices):
        BLANK = 0, _('Not annotated')
        IN_WORK = 1, _('In work')
        COMPLETE = 2, _('Validated')

    id = models.CharField(max_length=200, primary_key=True)
    species = models.CharField(max_length=200)
    strain = models.CharField(max_length=200)
    substrain = models.CharField(max_length=200)
    status = models.IntegerField(choices=Status.choices, default=Status.BLANK)
    
    def __str__(self):
        return self.id

class Chromosome (models.Model):
    
    id = models.CharField(max_length=200, primary_key=True)
    chromName = models.CharField(max_length=200)
    startPos = models.IntegerField(default=1)
    endPos = models.IntegerField()
    
    idGenome = models.ForeignKey(Genome, on_delete=models.CASCADE)
    
    class Meta:
        constraints = [
            models.UniqueConstraint("idGenome", "chromName", name="unique_name_of_chromosome_in_genome")
        ]

class ChromosomeSeq (models.Model):
    
    id = models.CharField(max_length=200, primary_key=True)
    sequence = models.CharField(max_length=1e10)
    
    idChrom = models.ForeignKey(Chromosome, on_delete=models.CASCADE)

class User (models.Model):
    
    class Role(models.IntegerChoices):
        READER = 0, _('Reader')
        ANNOTATOR = 1, _('Annotator')
        VALIDATOR = 2, _('Validator')
        
    email = models.CharField(max_length=200, primary_key=True)
    firstName = models.CharField(max_length=50)
    lastName = models.CharField(max_length=50)
    password = models.CharField(max_length=15)
    phoneNumber = models.IntegerField(max_length=12)
    role = models.IntegerField(choices=Role.choices, default=Role.READER)
    lastConnexion = models.DateField()

class Gene (models.Model):
    
    class Strand(models.IntegerChoices):
        SENSE = 1, _('Sense')
        ANTISENSE = -1, _('Antisense')
        
    class Status(models.IntegerChoices):
        ASSIGNABLE = 0, _('Assignable')
        BEING_ANNOTATED = 1, _('Being annotated')
        BEING_CORRECTED = 2, _('Being corrected')
        SUBMITTED = 3, _('Submitted to a validator')
        VALIDATED = 4, _('Validated')
    
    id = models.CharField(max_length=200, primary_key=True)
    geneName = models.CharField(max_length=200)
    geneSymbol = models.CharField(max_length=100)
    geneBiotype = models.CharField(max_length=200)
    strand = models.IntegerField(choices=Strand.choices, default=Strand.SENSE)
    startPos = models.IntegerField()
    endPos = models.IntegerField()
    description = models.CharField(max_length=700)
    status = models.IntegerField(choices=Status.choices, default=Status.ASSIGNABLE)
    
    idChrom = models.ForeignKey(Chromosome, on_delete=models.CASCADE)
    emailAnnotator = models.ForeignKey(User, on_delete=models.CASCADE)
    emailValidator = models.ForeignKey(User, on_delete=models.CASCADE)
    
class NucleotidicSeq (models.Model):
    id = models.CharField(max_length=200, primary_key=True)
    sequence = models.CharField(max_length=1e8) # au max : 99 999 999 pb
    
    idGene = models.ForeignKey(Gene, on_delete=models.CASCADE)
    
class Message(models.Model):
    
    # id auto-généré par django
    text = models.CharField(max_length=500)
    date = models.DateField(auto_now_add=True) # Automatically set the field to now when the object is first created
    
    emailAuthor = models.ForeignKey(User, on_delete=models.CASCADE)
    idGene = models.ForeignKey(Gene, on_delete=models.CASCADE)

class Peptide(models.Model):
    id = models.CharField(max_length=200, primary_key=True)
    transcriptName = models.CharField(max_length=200)
    transcriptBiotype = models.CharField(max_length=200)
    
    idGene = models.ForeignKey(Gene, on_delete=models.CASCADE)

class PeptideSeq(models.Model):
    id = models.CharField(max_length=200, primary_key=True)
    sequence = models.CharField(max_length=1e6) # au max : 999 999 acides aminés
    
    idPeptide = models.ForeignKey(Peptide, on_delete=models.CASCADE)

# https://docs.djangoproject.com/en/5.0/ref/models/fields/
# https://docs.djangoproject.com/en/5.0/topics/db/models/