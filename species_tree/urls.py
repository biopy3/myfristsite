from django.urls import path,re_path
from . import views

urlpatterns = [
    r'^plist/$'
    path(r'^result/download/$',views.download_results),
    path('manual.pdf',views.document),
    path('', views.home_page),
    path('save_post/',views.save_post),

]
