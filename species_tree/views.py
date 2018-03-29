from django.shortcuts import render,render_to_response
from .models import Records
from django.http import HttpResponse,HttpResponseRedirect,FileResponse
#from django.template import loader
from datetime import datetime
from .tasks import generate_tree
from django import forms
import os,string,random,uuid

class UserInfo(forms.Form):
    modelist = [['K80','K80'],['JC69','JC69'],['TS','TS'],['TV','TV'],['N','N'],['raw','raw'],
                ['F81','F81'],['K81','K81'],['F84','F84'],['BH87','BH87'],['T92','T92'],['TN93','TN93']
                ,['GG95','GG95'],['logdet','logdet'],['paralin','paralin'],['indel','indel'],
                ['indelblock','indelblock']]
    user = forms.CharField(widget=forms.TextInput(attrs={'class': 'special','size': '15'}),
                           max_length=16,min_length=3,error_messages={'required':u'user cannot be empty '})
    email = forms.EmailField(widget=forms.EmailInput(attrs={'class': 'special','size': '25'}),
                             error_messages={'required':u'please full in correct email'})
    inputfile = forms.FileField(widget=forms.FileInput(attrs={'style':'color:red'}),
                                error_messages={'required':u'please full in correct file'})
    model = forms.ChoiceField(label='model',choices=modelist)


modelist = ["K80","JC69","TS","TV","N","raw","F81","K81","F84","BH87", "T92","TN93","GG95","logdet","paralin","indel","indelblock"]
userinfo = UserInfo()

def document(request):
    cwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    the_file_name = 'manual.pdf'  # 显示在弹出对话框中的默认的下载文件名
    fname = '/species_tree/static/document/Program_of_identificate_species_dividing_line_manual.pdf'
    response = FileResponse(open(cwd +fname,'rb'))
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format(the_file_name)
    return response

def home_page(request):
    return render(request,"home.html",{'userinfo': userinfo,'modelist': modelist})

def download_results(request,access_code):
    try:
        record_ = Records.objects.get(access_code=access_code)
        fname = record_.resultfile.path
        the_file_name = fname.split('/')[-1]  # 显示在弹出对话框中的默认的下载文件名
        response = FileResponse(open(fname,'rb'))
        response['Content-Type'] = 'application/octet-stream'
        response['Content-Disposition'] = 'attachment;filename="{0}"'.format(the_file_name)
        return response
    except:
        return HttpResponse("Some errors hanppened!")

def save_post(request):
    success_str = "Submit successfully,waiting for minites we will send results to your email!"
    if request.method == "POST":
        user_input = UserInfo(request.POST,request.FILES) #request.POST is include all data
        if user_input.is_valid():
            user_name = user_input.cleaned_data['user']
            email = user_input.cleaned_data['email']
            inputfile = request.FILES['inputfile']
            model = user_input.cleaned_data['model']
            access_code = str(uuid.uuid1())
            record = Records.objects.create(user=user_name, inputfile=inputfile,
                                            access_code=access_code,email=email)
            infile_path = record.inputfile.path
            generate_tree.delay(infile_path,email,user_name,access_code,model)
            return HttpResponse("Thank you,we will send email for you!")
        else :
            error_msg = user_input.errors
            return render(request,'home.html',{'userinfo':userinfo,'errors':error_msg})
    else:
        return HttpResponse("errors")
