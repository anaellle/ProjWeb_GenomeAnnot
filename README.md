# ProjWeb_GenomeAnnot
__Emeline Bruyère<sup>1</sup>, Anaëlle Cossard<sup>1</sup>, Elora Vigo<sup>1</sup> et Alexis Michalowski-Skarbek<sup>1</sup>__
<br>
<sub>1. Université Paris-Saclay

Master 2 Web Project for functional annotation and analysis of bacterial genomes

### Requirements
Django version 4.2.7. <br>

- biopython==1.83
  
- django-bootstrap3==23.6
  
- django-filter==23.5
  
- django-tables2==2.7.0

- plotly==5.18.0 

-> See the requiements.txt file for details

### How-to use


- First be sure to clean the environment in ProjWeb_GenomeAnnot:
```
git clean -n -d -x
``` 

```
git clean -f -d -x
```
  
- Go to the GenomeAnnot directory :
```
cd GenomeAnnot
```
  - To create the database :
```
python manage.py makemigrations main 
``` 
then
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
