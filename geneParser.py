import random, math, string, sys, os
import matplotlib.pyplot as plt
from Bio import SeqIO

def CDSParser(filename):
    ### CREATION DU DICTIONNAIRE A FAIRE

    for record in SeqIO.parse("example.fasta", "fasta"):
        #### REMPLISSAGE DU DICTIONNAIRE
        pass
    
    return 0

def PepParser(filename):
    return 0

def GenomeParser(filename):
    return 0


def main(file):
    # On autorise aussi les .fasta j'imagine A FAIRE
    if(file.endswith('cds.fa')):
        CDSParser(file)
    elif(file.endswith('pep.fa')):
        PepParser(file)
    elif(file.endswith('.fa')):
        GenomeParser(file)
    else:
        # Est-ce qu'on gère les erreurs ici ? 
        print("file is not a fasta file")
    
    return 0