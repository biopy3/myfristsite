from django.urls import path,re_path
from . import views

urlpatterns = [
    path('result/download/',view.get_results)
    path('result/',views.query_get_results),
    path('manual.pdf',views.document),
    path('', views.home_page),
    path('save_post/',views.save_post),

]
