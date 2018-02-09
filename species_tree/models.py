from django.db import models

# Create your models here.

class Records(models.Model):
    user = models.CharField(max_length=200)
    inputfile = models.FileField(upload_to='recordsfile/')
    alignmentfile = models.FileField(upload_to='recordsfile/')
    phyfile = models.FileField(upload_to='recordsfile/')
    bar_chartfile = models.FileField(upload_to='recordsfile/')
    nwkfile = models.FileField(upload_to='recordsfile/')
    outputfile = models.FileField(upload_to='recordsfile/')
    submit_date = models.DateTimeField()
    email = models.EmailField()