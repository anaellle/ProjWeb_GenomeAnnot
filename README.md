# ProjWeb_GenomeAnnot
__Emeline Bruyère<sup>1</sup>, Anaëlle Cossard<sup>1</sup>, Elora Vigo<sup>1</sup> and Alexis Michalowski-Skarbek<sup>1</sup>__
<br>
<sub>1. Université Paris-Saclay

Master 2 Web Project for functional annotation and analysis of bacterial genomes. The goal was to produce a website allowing users to create an account with different roles and differents actions associted. The website enables to look at the database, visualize genomes, get the list of genes and peptides associated etc. Users can upload information to the database using a formular with fasta files.

### Requirements
Django version 4.2.7. <br>

- biopython==1.83
  
- django-bootstrap3==23.6
  
- django-filter==23.5
  
- django-tables2==2.7.0

- plotly==5.18.0 

-> See the requiements.txt file for details

### How-to use

- Go to the GenomeAnnot directory :
```
cd GenomeAnnot
```
  - To create the database :
```
python manage.py migrate
```
  - To init the database with 2 genomes :
```
python manage.py init-my-data
```
  - To run the website :
```
python manage.py runserver
```
