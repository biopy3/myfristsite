from django.shortcuts import render,render_to_response
from .models import Records
from django.http import HttpResponse,HttpResponseRedirect
#from django.template import loader
from datetime import datetime
from .tasks import generate_tree,modifytree
from django import forms

class UserInfo(forms.Form):
    user = forms.CharField(error_messages={'required':u'user cannot be empty '})
    email = forms.EmailField(error_messages={'required':u'please full in correct email'})
    inputfile = forms.FileField(error_messages={'required':u'please full in correct file'})

userinfo = UserInfo()


def display(request):
    return render(request,"display.html")
def home_page(request):
    return render(request,"home.html",{'userinfo': userinfo})

def save_post(request):
    success_str = "Submit successfully,waiting for minites we will send results to your email!"
    if request.method == "POST":
        user_input = UserInfo(request.POST,request.FILES) #request.POST is include all data
        if user_input.is_valid():
            user_name = user_input.cleaned_data['user']
            email = user_input.cleaned_data['email']
            inputfile = request.FILES['inputfile']
            record = Records.objects.create(user=user_name, inputfile=inputfile,
                                            submit_date=datetime.now(), email=email)
            infile_path = record.inputfile.path
            file_name = record.inputfile.name.split('.')[0]
            #generate_tree.delay(file_name,infile_path)
            modifytree.delay(file_name,infile_path,email,user_name)
            return HttpResponse(success_str)
        else :
            error_msg = user_input.errors
            return render(request,'home.html',{'userinfo':user_input,'errors':error_msg})

    return HttpResponse("errors")