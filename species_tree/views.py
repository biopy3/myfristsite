from django.shortcuts import render,render_to_response
from .models import Records
from django.http import HttpResponse,HttpResponseRedirect,FileResponse
#from django.template import loader
from datetime import datetime
from .tasks import generate_tree
from django import forms
import os,string,random

class UserInfo(forms.Form):
    user = forms.CharField(widget=forms.TextInput(attrs={'class': 'special','size': '15'}),
                           max_length=16,min_length=3,error_messages={'required':u'user cannot be empty '})
    email = forms.EmailField(widget=forms.EmailInput(attrs={'class': 'special','size': '25'}),
                             error_messages={'required':u'please full in correct email'})
    inputfile = forms.FileField(widget=forms.FileInput(attrs={'style':'color:red'}),
                                error_messages={'required':u'please full in correct file'})

userinfo = UserInfo()

class ResultInfo(forms.Form):
    email = forms.EmailField(widget=forms.TextInput(attrs={'class': 'special','size': '15'}),
                           error_messages={'required':u'user cannot be empty'})
    access_code = forms.CharField(widget=forms.TextInput(attrs={'style':'color:red'}),
                                  max_length=15,min_length=15,error_messages={'required':u'please full in correct file'})

resultinfo = ResultInfo()

def document(request):
    cwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    the_file_name = 'manual.pdf'  # 显示在弹出对话框中的默认的下载文件名
    fname = '/species_tree/static/document/Program_of_identificate_species_dividing_line_manual.pdf'
    response = FileResponse(open(cwd +fname,'rb'))
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format(the_file_name)
    return response

def home_page(request):
    return render(request,"home.html",{'userinfo': userinfo})

def query_get_results(request):
    if request.method == "post":
        query_info =  ResultInfo(request.POST)
        if query_info.is_valid():
            email = query_info.cleaned_data['email']
            access_code = query_info.cleaned_data['access_code']
            return HttpResponse(access_code)
            try:
                record_ = Records.objects.get(access_code=access_code)
                if email == record_.email:
                    fname = record_.resultfile.path
                    the_file_name = fname.split('/')[-1]  # 显示在弹出对话框中的默认的下载文件名
                    response = FileResponse(open(fname,'rb'))
                    response['Content-Type'] = 'application/octet-stream'
                    response['Content-Disposition'] = 'attachment;filename="{0}"'.format(the_file_name)
                    return response
                else:
                    return HttpResponse("Please put in correct infomation!")
            except:
                return HttpResponse("Please put in correct infomation!")
        else:
            error_msg = query_info.errors
            return render(request,'dislay.html',{'resultinfo':resultinfo,'errors':error_msg})

def save_post(request):
    success_str = "Submit successfully,waiting for minites we will send results to your email!"
    if request.method == "POST":
        user_input = UserInfo(request.POST,request.FILES) #request.POST is include all data
        if user_input.is_valid():
            user_name = user_input.cleaned_data['user']
            email = user_input.cleaned_data['email']
            inputfile = request.FILES['inputfile']
            access_code = ''.join(random.choice(string.digits + string.ascii_letters + string.punctuation) for _ in range(15))
            record = Records.objects.create(user=user_name, inputfile=inputfile,
                                            access_code=access_code,email=email)
            infile_path = record.inputfile.path
            generate_tree.delay(infile_path,email,user_name,access_code)
            return render(request,'display.html',{'resultinfo':resultinfo})
        else :
            error_msg = user_input.errors
            return render(request,'home.html',{'userinfo':resultinfo,'errors':error_msg})
    else:
        return HttpResponse("errors")
