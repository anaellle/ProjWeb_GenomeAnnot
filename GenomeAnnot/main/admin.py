from django.contrib import admin

# Register your models here.

from .models import *

admin.site.register(Genome)
admin.site.register(Chromosome)
admin.site.register(ChromosomeSeq)
admin.site.register(User)
admin.site.register(Gene)
admin.site.register(NucleotidicSeq)
admin.site.register(Message)
admin.site.register(Peptide)
admin.site.register(PeptideSeq)
