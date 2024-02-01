from fileParser import file_to_dico

import os
os.environ.setdefault('DJANGO_SETTINGS_MODULE','GenomeAnnot.settings')

from .models import Genome, Chromosome, ChromosomeSeq, Gene, NucleotidicSeq, Peptide, PeptideSeq

def fillDB(dico, table):
    for i in range(len(dico.keys())):
        table.objects.create(**dico[i])
        print('line : ',i,' added')

def addData(genomeDict, geneDict, pepDict):
    print('Starting genome dict')
    fillDB(genomeDict["genome"], Genome)
    print('Start of chromosome')
    fillDB(genomeDict["chromosome"], Chromosome)
    print('Start of chromosome sequence')
    fillDB(genomeDict["sequence"], ChromosomeSeq)
    print('start of gene')
    fillDB(geneDict["gene"], Gene)
    print('start of gene sequence')
    fillDB(geneDict["sequence"], NucleotidicSeq)
    print('start of peptide')
    fillDB(pepDict["peptide"], Peptide)
    print('start of peptide sequence')
    fillDB(pepDict["sequence"], PeptideSeq)


def test1():
    genome = file_to_dico('Escherichia_coli_str_k_12_substr_mg1655.fa')
    gene = file_to_dico('Escherichia_coli_str_k_12_substr_mg1655_cds.fa')
    peptide = file_to_dico('Escherichia_coli_str_k_12_substr_mg1655_pep.fa')
    addData(genome, gene, peptide)