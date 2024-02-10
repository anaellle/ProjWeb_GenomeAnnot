from fileParser import file_to_dico

import os

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "GenomeAnnot.settings")

from .models import Genome, Chromosome, ChromosomeSeq, Gene, NucleotidicSeq, Peptide, PeptideSeq


def fillDBGenome(dico):
    for id in list(dico.keys()):
        Genome.objects.create(**dico[id])


def fillDBChromosome(dico):
    for id in list(dico.keys()):
        current_line = dico[id]
        genomeID = current_line.pop("idGenome")
        genomeRef = Genome.objects.get(pk=genomeID)
        Chromosome.objects.create(**current_line, idGenome=genomeRef)


def fillDBChromosomeSeq(dico):
    for id in list(dico.keys()):
        current_line = dico[id]
        chromosomeID = current_line.pop("idChrom")
        chromRef = Chromosome.objects.get(pk=chromosomeID)
        ChromosomeSeq.objects.create(**current_line, idChrom=chromRef)


def fillDBGene(dico):
    for id in list(dico.keys()):
        current_line = dico[id]
        chromosomeID = current_line.pop("idChrom")
        chromRef = Chromosome.objects.get(pk=chromosomeID)
        Gene.objects.create(**current_line, idChrom=chromRef)


def fillDBGeneSeq(dico):
    for id in list(dico.keys()):
        current_line = dico[id]
        geneID = current_line.pop("idGene")
        geneRef = Gene.objects.get(pk=geneID)
        NucleotidicSeq.objects.create(**current_line, idGene=geneRef)


def fillDBPep(dico):
    for id in list(dico.keys()):
        current_line = dico[id]
        geneID = current_line.pop("idGene")
        geneRef = Gene.objects.get(pk=geneID)
        Peptide.objects.create(**current_line, idGene=geneRef)


def fillDBPepSeq(dico):
    for id in list(dico.keys()):
        current_line = dico[id]
        pepID = current_line.pop("idPeptide")
        pepRef = Peptide.objects.get(pk=pepID)
        PeptideSeq.objects.create(**current_line, idPeptide=pepRef)


def addData(genomeDict, geneDict, pepDict):
    ## Check annotations
    genome = list(genomeDict["genome"].keys())[0]
    if geneDict['status'] == len(geneDict["gene"].keys()):
        genomeDict["genome"][genome]["status"]=2
    elif 0 < geneDict['status'] < len(geneDict["gene"].keys()):
        genomeDict["genome"][genome]["status"]=1

    ## Ajout du génome
    print("Adding genome infos")
    fillDBGenome(genomeDict["genome"])

    ## Ajout des chromosomes
    print("Adding chromosomes infos")
    fillDBChromosome(genomeDict["chromosome"])

    ## Ajout des sequences des chromosomes
    print("Adding chromosomes sequences")
    fillDBChromosomeSeq(genomeDict["sequence"])

    ## Ajout des gènes
    print("Adding genes")
    fillDBGene(geneDict["gene"])

    ## Ajout des séquences nucléotidiques
    print("Adding nucleotidic sequences")
    fillDBGeneSeq(geneDict["sequence"])

    ## Ajout des peptides
    print("Adding peptides")
    fillDBPep(pepDict["peptide"])

    ## Ajout des séquences peptidiques
    print("Adding peptides sequences")
    fillDBPepSeq(pepDict["sequence"])


def downloadAndFill(genomeFile, geneFile, peptideFile):
    genome = file_to_dico(genomeFile)
    gene = file_to_dico(geneFile)
    peptide = file_to_dico(peptideFile)
    addData(genome, gene, peptide)
