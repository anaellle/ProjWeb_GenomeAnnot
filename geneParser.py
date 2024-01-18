import sys, os
from Bio import SeqIO

def CDSParser(filename):
    # création du dictionnaire
    CDSdict = {}

    # on boucle sur toutes les séquences
    for record in SeqIO.parse(filename, "fasta"):
        # création du dictionnaire pour la séquence
        CDSdict[record.id] = {}
        # ajout de la séquence
        CDSdict[record.id]['sequence'] = str(record.seq)

        infos = record.description.split()
        chr = infos[2].split(':')

        # ajout des autres informations
        CDSdict[record.id]['geneName'] = infos[3].split(':')[1]
        CDSdict[record.id]['geneSymbol'] = infos[6].split(':')[1]
        CDSdict[record.id]['geneBiotype'] = infos[4].split(':')[1]
        CDSdict[record.id]['chromosome'] = chr[1]
        CDSdict[record.id]['strand'] = chr[5]
        CDSdict[record.id]['startPos'] = chr[3]
        CDSdict[record.id]['endPos'] = chr[4]
        CDSdict[record.id]['description'] =  record.description.split('description:')[1]
    
    return CDSdict

def PepParser(filename):
    pepdict={}

    for record in SeqIO.parse(filename, "fasta"):
        pepdict[record.id] = {}
        pepdict[record.id]['sequence'] = str(record.seq)

        infos = record.description.split()

        pepdict[record.id]['transcritName'] = infos[4].split(':')[1]
        pepdict[record.id]['transcritBiotype'] = infos[6].split(':')[1]
        pepdict[record.id]['description'] =  record.description.split('description:')[1]
    
    return pepdict

def GenomeParser(filename):
    genome = list(SeqIO.parse(filename, "fasta"))

    gendict = {}

    gendict["genomeID"] = genome[0].description.split()[2].split(':')[1]
    gendict["species"] = filename.split('.')[0]
    gendict["strain"] = genome[0].description.split()[2].split(':')[5]
    gendict["sequence"] = str(genome[0].seq)
    
    return gendict


def main(file):
    # On autorise aussi les .fasta ??q
    if(file.endswith('cds.fa')):
        res = CDSParser(file)
    elif(file.endswith('pep.fa')):
        res = PepParser(file)
    elif(file.endswith('.fa')):
        res = GenomeParser(file)
    else:
        # Est-ce qu'on gère les erreurs ici ? 
        print("file is not a fasta file")
    
    return res