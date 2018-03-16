from django.db import models

# Create your models here.

class Records(models.Model):
    user = models.CharField(max_length=200)
    inputfile = models.FileField(upload_to ='species_tree/recordsfile/')
    access_code = models.CharField(max_length=15,min_length=15)
    resultfile = models.FileField()
    submit_date = models.DateTimeField()
    email = models.EmailField()
