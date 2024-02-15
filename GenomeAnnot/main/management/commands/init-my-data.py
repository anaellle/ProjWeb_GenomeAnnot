from django.core.management.base import BaseCommand
from main.insertion import uploadAndFill
import os

class Command(BaseCommand):
    help = 'Initial filling of the database with files in the Data repository'
    
    def handle(self, *args, **options):
        '''
        Fill the database with the fasta files available in the 
        directory GenomeAnnot/Data
        '''
        current_dir = os.getcwd()
        data_dir = os.path.join(current_dir, 'Data')

        if os.path.exists(data_dir) and os.path.isdir(data_dir):
            # Get the list of files in the directory
            files_in_directory = os.listdir(data_dir)
            
            # loop through files to get the genome, gene and peptide files together
            for file in files_in_directory:
                name = file.split('.fa')[0]
                
                # get the name of the other files
                if name.endswith('_cds'):
                    genomeName = name.split('_cds')[0]+'.fa'
                    geneName = name+'.fa'
                    peptideName = name.split('cds')[0]+'pep.fa'
                
                elif name.endswith('_pep'):
                    genomeName = name.split('_pep')[0]+'.fa'
                    geneName = name.split('pep')[0]+'cds.fa'
                    peptideName = name+'.fa'
                
                else :
                    genomeName = name+'.fa'
                    geneName = name+'_cds.fa'
                    peptideName = name+'_pep.fa'
                
                # get the paths of the files
                genomePath = os.path.join(data_dir,genomeName)
                genePath = os.path.join(data_dir,geneName)
                peptidePath = os.path.join(data_dir,peptideName)
                
                # check if they exists and add them to the database
                if os.path.exists(genomePath) and os.path.exists(genePath) and os.path.exists(peptidePath):
                    self.stdout.write(self.style.WARNING('Adding {} genome'.format(genomeName)))
                    uploadAndFill(genomePath,genePath,peptidePath)
                    files_in_directory.remove(genomeName)
                    files_in_directory.remove(geneName)
                    files_in_directory.remove(peptideName)

            self.stdout.write(self.style.SUCCESS('Your files have been successfully added to the database'))
        else :
            self.stdout.write(self.style.ERROR('Nothing as been added to the database'))


