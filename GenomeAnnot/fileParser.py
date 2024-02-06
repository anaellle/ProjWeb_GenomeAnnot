from Bio import SeqIO

def CDSParser(filename):
    # création des dictionnaires
    CDSdict = {}
    CDSdict["gene"] = {}
    CDSdict["sequence"] = {}

    # on boucle sur toutes les séquences
    for record in SeqIO.parse(filename, "fasta"):
        # création des dictionnaires pour la séquence et le gène
        CDSdict["gene"][record.id] = {}
        CDSdict["sequence"][record.id] = {}
        
        ## ajout des informations de la séquence
        CDSdict["sequence"][record.id]['id'] = record.id
        CDSdict["sequence"][record.id]['idGene'] = record.id
        CDSdict["sequence"][record.id]['sequence'] = str(record.seq)

        ## ajout des informations du gène
        # récupération des informations
        infos = record.description.split()
        attributes = ['gene', 'gene_biotype', 'gene_symbol']
        names = ['geneName', 'geneBiotype', 'geneSymbol']

        CDSdict["gene"][record.id]['id'] = record.id 

        # stockage des informations du gène
        for att in range(len(attributes)) :
            index = [i for i in range(len(infos)) if attributes[att] in infos[i]]
            if index != []:
                CDSdict["gene"][record.id][names[att]] = infos[index[0]].split(':')[1]
            elif attributes[att] == 'gene':
                CDSdict["gene"][record.id][names[att]] = record.id

        # information du chromosome
        index_pos = [i for i in range(len(infos)) if 'chromosome' in infos[i]]
        chr_infos = infos[index_pos[0]].split(':')
        CDSdict["gene"][record.id]['idChrom'] = chr_infos[1]
        CDSdict["gene"][record.id]['startPos'] = chr_infos[3]
        CDSdict["gene"][record.id]['endPos'] = chr_infos[4]
        if len(chr_infos) > 5 :
            CDSdict["gene"][record.id]['strand'] = chr_infos[5]

        # ajout de la description
        if 'description' in record.description:
            CDSdict["gene"][record.id]['descriptionGene'] =  record.description.split('description:')[1]
    
    return CDSdict


def PepParser(filename):
    pepdict={}
    pepdict["peptide"] = {}
    pepdict["sequence"] = {}

    for record in SeqIO.parse(filename, "fasta"):
        pepdict["peptide"][record.id] = {}
        pepdict["sequence"][record.id] = {}
        
        # dictionnaire pour la séquence
        pepdict["sequence"][record.id]['id'] = record.id
        pepdict["sequence"][record.id]['idPeptide'] = record.id
        pepdict["sequence"][record.id]['sequence'] = str(record.seq)

        # dictionnaire pour le peptide
        infos = record.description.split()
        attributes = ['transcript', 'transcript_biotype']
        names = ['transcriptName', 'transcriptBiotype']

        pepdict["peptide"][record.id]['id'] = record.id
        pepdict["peptide"][record.id]['idGene'] = record.id

        # ajout des informations du peptide
        for att in range(len(attributes)) :
            index = [i for i in range(len(infos)) if attributes[att] in infos[i]]
            if index != []:
                pepdict["peptide"][record.id][names[att]] = infos[index[0]].split(':')[1]
            elif attributes[att] == 'transcript' :
                pepdict["peptide"][record.id][names[att]] = record.id

        # ajout de la description
        if 'description' in record.description:
            pepdict["peptide"][record.id]['descriptionPep'] =  record.description.split('description:')[1]
    
    return pepdict


def GenomeParser(filename):
    # création des dictionnaires
    gendict = {}
    gendict["genome"] = {}
    gendict["chromosome"] = {}
    gendict["sequence"] = {}

    # Nom du fichier
    genomeID = filename.split('.')[0]

    # Formatage du nom
    strain = genomeID.find('_str')
    substr = genomeID.find('_substr')

    flagSTR = True
    flagSUB = True

    if strain != -1:
        species = genomeID[:strain]
        if substr != -1 :
            strain = genomeID[strain+5:substr]
            substrain = genomeID[substr+6:]
            genomeName = species+'_STR_'+strain+'_SUBSTR_'+substrain
        else :
            strain = genomeID[strain+5:]
            genomeName = species+'_STR_'+strain
            flagSUB = False
    else :
        species = genomeID
        genomeName = genomeID
        flagSTR = False
        flagSUB = False

    # Ajout des informations du génome
    gendict["genome"][genomeName]={}
    gendict["genome"][genomeName]['id'] = genomeName
    gendict["genome"][genomeName]['species'] = species

    if flagSTR:
        gendict["genome"][genomeName]['strain'] = strain
    if flagSUB:
        gendict["genome"][genomeName]['substrain'] = substrain

    # Remplissage pour les chromosomes et les sequences
    for record in SeqIO.parse(filename, "fasta"):
        chr_infos = record.description.split()[2]
        chrName = chr_infos.split(':')[1]

        # Ajout de la sequence
        gendict["sequence"][chrName] = {}
        gendict["sequence"][chrName]['id'] = chrName
        gendict["sequence"][chrName]['idChrom'] = chrName
        gendict["sequence"][chrName]['sequence'] = str(record.seq)

        # Ajout des informations du chromosome
        gendict["chromosome"][chrName] = {}
        
        gendict["chromosome"][chrName]['id'] = chrName
        gendict["chromosome"][chrName]['chromName'] = chrName
        gendict["chromosome"][chrName]["startPos"] = chr_infos.split(':')[3]
        gendict["chromosome"][chrName]["endPos"] = chr_infos.split(':')[4]
        if len(chr_infos) > 5 :
            gendict["chromosome"][chrName]["strand"] = record.description.split()[2].split(':')[5]
        
        gendict["chromosome"][chrName]['idGenome'] = genomeName

    return gendict


def file_to_dico(file):
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