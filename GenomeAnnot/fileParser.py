from Bio import SeqIO
from django.core.files import File

def CDSParser(filename):
    '''
    Convert a gene fasta file into a dictionary.
    It contains 3 keys :
        - gene : a dictionary, for each gene another dictionary is associated
        with his available attributes
        - sequence : a dictionary, for each gene his sequence is put in a dictionnary
        - status : number of gene considered as annotated
    '''
    # dictionaries creation
    CDSdict = {}
    CDSdict["gene"] = {}
    CDSdict["sequence"] = {}

    # to know how many genes are annotated
    cpt = 0

    # check the format of the file
    if isinstance(filename, File):
        records = SeqIO.parse(filename.temporary_file_path(), "fasta")
    elif isinstance(filename, str):
        records = SeqIO.parse(filename, "fasta")

    # loop trhough the genes
    for record in records:
        # dictionaries creation for the gene and the sequence
        CDSdict["gene"][record.id] = {}
        CDSdict["sequence"][record.id] = {}
        
        ## add sequence information
        CDSdict["sequence"][record.id]['id'] = record.id
        CDSdict["sequence"][record.id]['idGene'] = record.id
        CDSdict["sequence"][record.id]['sequence'] = str(record.seq)

        ## add gene information
        # get the info
        infos = record.description.split()
        attributes = ['gene', 'gene_biotype', 'gene_symbol']
        names = ['geneName', 'geneBiotype', 'geneSymbol']

        CDSdict["gene"][record.id]['id'] = record.id 

        # store info if available
        for att in range(len(attributes)) :
            index = [i for i in range(len(infos)) if attributes[att] in infos[i]]
            if index != []:
                CDSdict["gene"][record.id][names[att]] = infos[index[0]].split(':')[1]
            elif attributes[att] == 'gene':
                CDSdict["gene"][record.id][names[att]] = record.id

        # chromosome information
        index_pos = [i for i in range(len(infos)) if ('chromosome' in infos[i] or 'plasmid' in infos[i])]
        chr_infos = infos[index_pos[0]].split(':')
        CDSdict["gene"][record.id]['idChrom'] = chr_infos[1]
        CDSdict["gene"][record.id]['startPos'] = chr_infos[3]
        CDSdict["gene"][record.id]['endPos'] = chr_infos[4]
        if len(chr_infos) > 5 :
            CDSdict["gene"][record.id]['strand'] = chr_infos[5]

        # add description if available and consider gene as annotated
        if 'description' in record.description:
            CDSdict["gene"][record.id]['descriptionGene'] =  record.description.split('description:')[1]
            CDSdict["gene"][record.id]['status'] = 4
            cpt = cpt +1

        CDSdict["status"] = cpt
    
    return CDSdict


def PepParser(filename):
    '''
    Convert a peptide fasta file into a dictionary.
    It contains 2 keys :
        - peptide : a dictionary, for each peptide another dictionary is associated
        with his available attributes
        - sequence : a dictionnary, for each peptide his sequence is put in a dictionary
    '''
    # dictionaries creation
    pepdict={}
    pepdict["peptide"] = {}
    pepdict["sequence"] = {}

    # check file
    if isinstance(filename, File):
        records = SeqIO.parse(filename.temporary_file_path(), "fasta")
    elif isinstance(filename, str):
        records = SeqIO.parse(filename, "fasta")

    # loop through each peptide
    for record in records:
        pepdict["peptide"][record.id] = {}
        pepdict["sequence"][record.id] = {}
        
        # sequence information
        pepdict["sequence"][record.id]['id'] = record.id
        pepdict["sequence"][record.id]['idPeptide'] = record.id
        pepdict["sequence"][record.id]['sequence'] = str(record.seq)

        # peptide information
        infos = record.description.split()
        attributes = ['transcript', 'transcript_biotype']
        names = ['transcriptName', 'transcriptBiotype']

        pepdict["peptide"][record.id]['id'] = record.id
        pepdict["peptide"][record.id]['idGene'] = record.id

        # add information if available
        for att in range(len(attributes)) :
            index = [i for i in range(len(infos)) if attributes[att] in infos[i]]
            if index != []:
                pepdict["peptide"][record.id][names[att]] = infos[index[0]].split(':')[1]
            elif attributes[att] == 'transcript' :
                pepdict["peptide"][record.id][names[att]] = record.id

        # add description if available
        if 'description' in record.description:
            pepdict["peptide"][record.id]['descriptionPep'] =  record.description.split('description:')[1]
    
    return pepdict


def GenomeParser(filename):
    '''
    Convert a genome fasta file into a dictionnary.
    It contains 3 key :
        - genome : a dictionary, contains all the information of the genome
        - chromosome : a dictionary, for each chromosome contains his information
        - sequence : a dictionary, for each chromosome his sequence is stored in a dictionary
    '''
    # dictionaries creation
    gendict = {}
    gendict["genome"] = {}
    gendict["chromosome"] = {}
    gendict["sequence"] = {}

    # Check file
    if isinstance(filename, str):
        genomeID = filename.split('.')[0].split('/')[-1]
    elif isinstance(filename, File):
        genomeID = filename.name.split('.')[0]

    # Creation of genome's name
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

    # add genome information
    gendict["genome"][genomeName]={}
    gendict["genome"][genomeName]['id'] = genomeName
    gendict["genome"][genomeName]['species'] = species

    if flagSTR:
        gendict["genome"][genomeName]['strain'] = strain
    if flagSUB:
        gendict["genome"][genomeName]['substrain'] = substrain

    if isinstance(filename, File):
        records = SeqIO.parse(filename.temporary_file_path(), "fasta")
    elif isinstance(filename, str):
        records = SeqIO.parse(filename, "fasta")

    # loop through each chromosome
    for record in records:
        chr_infos = record.description.split()[2]
        chrName = chr_infos.split(':')[1]

        # add sequence information
        gendict["sequence"][chrName] = {}
        gendict["sequence"][chrName]['id'] = chrName
        gendict["sequence"][chrName]['idChrom'] = chrName
        gendict["sequence"][chrName]['sequence'] = str(record.seq)

        # add chromosome information
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
    '''
    Take a fasta file in input and return a dictionary.
    Check if the file is the path or the actual file and send it to 
    the right parser.
    '''
    if isinstance(file, str):
        if(file.endswith('cds.fa')):
            res = CDSParser(file)
        elif(file.endswith('pep.fa')):
            res = PepParser(file)
        elif(file.endswith('.fa')):
            res = GenomeParser(file)
    
    elif isinstance(file, File):
        if(file.name.endswith('cds.fa')):
            res = CDSParser(file)
        elif(file.name.endswith('pep.fa')):
            res = PepParser(file)
        elif(file.name.endswith('.fa')):
            res = GenomeParser(file)
    
    else: 
        print("The file does not have the required format")
        return
    return res