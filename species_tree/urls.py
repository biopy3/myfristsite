from django.urls import path,re_path
from . import views

urlpatterns = [
    path('home_page/', views.home_page),
    path('',views.save_post),

]