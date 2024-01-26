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

        # récupération des informations
        infos = record.description.split()
        attributes = ['gene', 'gene_biotype', 'gene_symbol']
        names = ['geneName', 'geneBiotype', 'geneSymbol']

        # ajout des informations du gène
        for att in range(len(attributes)) :
            index = [i for i in range(len(infos)) if attributes[att] in infos[i]]
            if index != []:
                CDSdict[record.id][names[att]] = infos[index[0]].split(':')[1]
            else :
                CDSdict[record.id][names[att]] = 'NULL'

        # information du chromosome
        index_pos = [i for i in range(len(infos)) if 'chromosome' in infos[i]]
        if len(index_pos) > 0:
            chr_infos = infos[index_pos[0]].split(':')
            CDSdict[record.id]['idChrom'] = chr_infos[1]
            CDSdict[record.id]['startPos'] = chr_infos[3]
            CDSdict[record.id]['endPos'] = chr_infos[4]
            if len(chr_infos) > 5 :
                CDSdict[record.id]['strand'] = chr_infos[5]
            else :
                CDSdict[record.id]['strand'] = 'NULL'
        else :
            CDSdict[record.id]['idChrom'] = 'NULL'
            CDSdict[record.id]['startPos'] = 'NULL'
            CDSdict[record.id]['endPos'] = 'NULL'
            CDSdict[record.id]['strand'] = 'NULL'

        # ajout de la description
        if 'description' in record.description:
            CDSdict[record.id]['description'] =  record.description.split('description:')[1]
        else :
            CDSdict[record.id]['description'] = 'NULL'
    
    return CDSdict


def PepParser(filename):
    pepdict={}

    pepdict["IDgenome"] = filename.split('_pep')[0]

    for record in SeqIO.parse(filename, "fasta"):
        pepdict[record.id] = {}
        pepdict[record.id]['sequence'] = str(record.seq)

        infos = record.description.split()
        attributes = ['transcript', 'transcript_biotype']
        names = ['transcritName', 'transcritBiotype']

        # ajout des informations du gène
        for att in range(len(attributes)) :
            index = [i for i in range(len(infos)) if attributes[att] in infos[i]]
            if index != []:
                pepdict[record.id][names[att]] = infos[index[0]].split(':')[1]
            else :
                pepdict[record.id][names[att]] = 'NULL'

        # ajout de la description
        if 'description' in record.description:
            pepdict[record.id]['description'] =  record.description.split('description:')[1]
        else :
            pepdict[record.id]['description'] = 'NULL'
    
    return pepdict


def GenomeParser(filename):
    gendict = {}

    genomeID = filename.split('.')[0]

    gendict[genomeID]={}
    gendict["IDgenome"] = genomeID


    strain = genomeID.find('_str')
    substr = genomeID.find('_substr')

    if strain != -1:
        gendict[genomeID]["species"] = genomeID[:strain]
        if substr != -1 :
            gendict[genomeID]["strain"] = genomeID[strain+5:substr]
            gendict[genomeID]["substrain"] = genomeID[substr+6:]
            #genome_name = gendict[genomeID]["species"]+'_STR_'+gendict[genomeID]["strain"]+'_SUBSTR_'+gendict[genomeID]["substrain"]
        else :
            gendict[genomeID]["strain"] = genomeID[strain+5:]
            #genome_name = gendict[genomeID]["species"]+'_STR_'+gendict[genomeID]["strain"]
    else :
        gendict[genomeID]["species"] = genomeID
        gendict[genomeID]["strain"] = 'NULL'
        gendict[genomeID]["substrain"] = 'NULL'

    #gendict["IDgenome"] = genome_name
    
    gendict[genomeID]["status"] = 0

    for record in SeqIO.parse(filename, "fasta"):
        chr_infos = record.description.split()[2]
        chrName = chr_infos.split(':')[1]
        gendict[chrName] = {}
        
        gendict[chrName]["startPos"] = chr_infos.split(':')[3]
        gendict[chrName]["endPos"] = chr_infos.split(':')[4]
        if len(chr_infos) > 5 :
            gendict[chrName]["strand"] = record.description.split()[2].split(':')[5]
        else :
            gendict[chrName]['strand'] = 'NULL'

        gendict[chrName]["sequence"] = str(record.seq)
    
    return gendict


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
