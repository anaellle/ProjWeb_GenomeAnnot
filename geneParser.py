import sys, os
from Bio import SeqIO

def CDSParser(filename):
    # création du dictionnaire
    CDSdict = {}
    CDSdict["IDgenome"] = filename.split('_cds')[0]

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
        CDSdict[record.id]['idChrom'] = chr[1]
        CDSdict[record.id]['strand'] = chr[5]
        CDSdict[record.id]['startPos'] = chr[3]
        CDSdict[record.id]['endPos'] = chr[4]
        CDSdict[record.id]['description'] =  record.description.split('description:')[1]
    
    return CDSdict


def PepParser(filename):
    pepdict={}

    pepdict["IDgenome"] = filename.split('_pep')[0]

    for record in SeqIO.parse(filename, "fasta"):
        pepdict[record.id] = {}
        pepdict[record.id]['sequence'] = str(record.seq)

        infos = record.description.split()

        pepdict[record.id]['transcritName'] = infos[4].split(':')[1]
        pepdict[record.id]['transcritBiotype'] = infos[6].split(':')[1]
        pepdict[record.id]['description'] =  record.description.split('description:')[1]
    
    return pepdict


def GenomeParser(filename):
    gendict = {}
    chrdict = {}

    genomeID = filename.split('.')[0]

    gendict[genomeID]={}
    chrdict["IDgenome"] = genomeID

    gendict[genomeID]["species"] = "NULL"
    gendict[genomeID]["strain"] = "NULL"
    gendict[genomeID]["substrain"] = "NULL"
    gendict[genomeID]["status"] = 0

    for record in SeqIO.parse(filename, "fasta"):
        chrName = record.description.split()[2].split(':')[1]
        chrdict[chrName] = {}
        chrdict[chrName]["startPos"] = record.description.split()[2].split(':')[3]
        chrdict[chrName]["endPos"] = record.description.split()[2].split(':')[4]
        chrdict[chrName]["strand"] = record.description.split()[2].split(':')[5]
        chrdict[chrName]["sequence"] = str(record.seq)
    
    res = {}
    res["genome"] = gendict
    res["chromosome"] = chrdict

    return res


def main(file):
    if(file.endswith('cds.fa')):
        res = CDSParser(file)
    elif(file.endswith('pep.fa')):
        res = PepParser(file)
    elif(file.endswith('.fa')):
        res = GenomeParser(file)
    else: 
        print("file is not a .fa file")
        res = -1
    
    return res

