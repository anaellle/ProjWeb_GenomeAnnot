from django.core.management.base import BaseCommand
from main.insertion import uploadAndFill
from django.core.files import File

class Command(BaseCommand):
    help = 'Initial filling of the database with files'

    def add_arguments(self, parser):
        parser.add_argument('genomePath', type=str, help='Path of the genome file')
        parser.add_argument('genePath', type=str, help='Path of the gene file')
        parser.add_argument('peptidePath', type=str, help='Path of the peptide file')

    def handle(self, *args, **options):
        # Your function goes here
        genomePath = options['genomePath']

        genePath = options['genePath']

        peptidePath = options['peptidePath']

        uploadAndFill(genomePath,genePath,peptidePath)
        self.stdout.write(self.style.SUCCESS('Your files have been successfully added to the database'))

