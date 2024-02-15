from fileParser import file_to_dico
from django.db import transaction

import os

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "GenomeAnnot.settings")

from .models import Genome, Chromosome, ChromosomeSeq, Gene, NucleotidicSeq, Peptide, PeptideSeq


def fillDBGenome(dico):
    '''
    Takes the genome dictionary as input and fill the table Genome
    '''
    with transaction.atomic():
        for id in list(dico.keys()):
            Genome.objects.create(**dico[id])


def fillDBChromosome(dico):
    '''
    Takes the chromosome dictionary as input and fill the table Chromosome
    '''
    with transaction.atomic():
        for id in list(dico.keys()):
            current_line = dico[id]
            genomeID = current_line.pop("idGenome")
            genomeRef = Genome.objects.get(pk=genomeID)
            Chromosome.objects.create(**current_line, idGenome=genomeRef)


def fillDBChromosomeSeq(dico):
    '''
    Takes the chromosome's sequence dictionary as input and fill the table ChromosomeSeq
    '''
    with transaction.atomic():
        for id in list(dico.keys()):
            current_line = dico[id]
            chromosomeID = current_line.pop("idChrom")
            chromRef = Chromosome.objects.get(pk=chromosomeID)
            ChromosomeSeq.objects.create(**current_line, idChrom=chromRef)


def fillDBGene(dico):
    '''
    Takes the gene dictionary as input and fill the table Gene
    '''
    with transaction.atomic():
        for id in list(dico.keys()):
            current_line = dico[id]
            chromosomeID = current_line.pop("idChrom")
            chromRef = Chromosome.objects.get(pk=chromosomeID)
            Gene.objects.create(**current_line, idChrom=chromRef)


def fillDBGeneSeq(dico):
    '''
    Takes the genes sequences dictionary as input and fill the table NucleotificSeq
    '''
    with transaction.atomic():
        for id in list(dico.keys()):
            current_line = dico[id]
            geneID = current_line.pop("idGene")
            geneRef = Gene.objects.get(pk=geneID)
            NucleotidicSeq.objects.create(**current_line, idGene=geneRef)


def fillDBPep(dico):
    '''
    Takes the peptide dictionary as input and fill the table Peptide
    '''
    with transaction.atomic():
        for id in list(dico.keys()):
            current_line = dico[id]
            geneID = current_line.pop("idGene")
            geneRef = Gene.objects.get(pk=geneID)
            Peptide.objects.create(**current_line, idGene=geneRef)


def fillDBPepSeq(dico):
    '''
    Takes the peptides sequences dictionary as input and fill the table PeptideSeq
    '''
    with transaction.atomic():
        for id in list(dico.keys()):
            current_line = dico[id]
            pepID = current_line.pop("idPeptide")
            pepRef = Peptide.objects.get(pk=pepID)
            PeptideSeq.objects.create(**current_line, idPeptide=pepRef)


def addData(genomeDict, geneDict, pepDict):
    '''
    Takes as input 3 dictionaries for the genome, genes and peptides.
    Call the right function for each to add information in the database.
    '''

    ## Check genes annotation
    genome = list(genomeDict["genome"].keys())[0]
    if geneDict['status'] == len(geneDict["gene"].keys()):
        genomeDict["genome"][genome]["status"]=2
    elif 0 < geneDict['status'] < len(geneDict["gene"].keys()):
        genomeDict["genome"][genome]["status"]=1

    ## add genome
    print("Adding genome infos")
    fillDBGenome(genomeDict["genome"])

    ## add chromosome
    print("Adding chromosomes infos")
    fillDBChromosome(genomeDict["chromosome"])

    ## add chromosome's sequence
    print("Adding chromosomes sequences")
    fillDBChromosomeSeq(genomeDict["sequence"])

    ## add genes
    print("Adding genes")
    fillDBGene(geneDict["gene"])

    ## add genes sequences
    print("Adding nucleotidic sequences")
    fillDBGeneSeq(geneDict["sequence"])

    ## add peptides
    print("Adding peptides")
    fillDBPep(pepDict["peptide"])

    ## add peptides sequences
    print("Adding peptides sequences")
    fillDBPepSeq(pepDict["sequence"])


def uploadAndFill(genomeFile, geneFile, peptideFile):
    '''
    Main function : takes as input the 3 files
    Call the file_to_dico funtcion to tranform each file in a dictionary
    Call the addData function with the dictionaries to fill the database
    '''
    genome = file_to_dico(genomeFile)
    gene = file_to_dico(geneFile)
    peptide = file_to_dico(peptideFile)
    addData(genome, gene, peptide)
