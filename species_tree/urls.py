from django.urls import path,re_path
from . import views

urlpatterns = [
    path('clustalx_save/',views.clutalx_save),
    path('result_download/<uuid:access_code>',views.clustalx_result_download),
    path('clustalx/',views.clustalx_page),
    path('result/download/<uuid:access_code>',views.download_results),
    path('manual.pdf',views.document),
    path('', views.home_page),
    path('save_post/',views.save_post),

]
