import argparse
import sys
import warnings
import os
import datetime 
from dateutil import tz
import pandas as pd

# Suppress FutureWarnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)


#call local path 
#local_path = os.getcwd()
#code_root = os.path.abspath(local_path+'/src')

#if code_root not in sys.path:
#    sys.path.insert(0, code_root)
    

#Paths for jupyter
current_working_directory = os.getcwd()
#local_path = os.path.dirname(os.path.realpath(__file__))
projDir = os.path.abspath(os.path.join(current_working_directory, '../..', ''))
dataDir = os.path.abspath(os.path.join(current_working_directory, '..', 'Data'))
dataDirI = dataDir+'/Input'
dataDirO = dataDir+'/Output'

#Paths for python 
current_working_directory2 = os.getcwd()+'/NIRSE'
projDir2 = os.path.abspath(os.path.join(current_working_directory2, '../..', ''))
dataDir2 = os.path.abspath(os.path.join(current_working_directory2, '..', 'Data'))
dataDirI2 = dataDir2+'/Input'
dataDirO2 = dataDir2+'/Output'

current_working_directory3 = os.getcwd()
projDir3 = os.path.abspath(os.path.join(current_working_directory3, '../..', ''))
dataDir3 = os.path.abspath(os.path.join(current_working_directory3, '../..', 'Data'))
dataDirI3 = dataDir3+'/Input'
dataDirO3 = dataDir3+'/Output'

#Input
egresos_48 = '/egresos_full_age_48.parquet'
tribu_priv = '/EstadoCargaIEEH__SEREMI_18122024.xlsx'
tribu_pub = '/EstadoCargaIEEH_SS_181224.xlsx' 
defunciones = '/Ficha Activo Información Defunciones.xlsx'
hospitalaria = '/Ficha Activo Información Egresos_Hospitalarios 2.0.xlsx'
nacimientos = '/Ficha Activo Información Nacimientos.xlsx'
egresos2024 ='/2024-12-16_Egresos.csv'
egresos = '/egresos.csv'
establecimientos = '/establecimientos-full.xlsx'
data_nirse = '/data.csv'
vigilancia = '/EGRESOS_VIGILANCIA_all_cause_under5_encr.csv'
nacimientos_2024 = '/NAC_RNI_EGRESOS_ENTREGA_ISCI_20_12_2024_encr.csv'
urgencias = [f'/DEIS/AtencionesUrgencia20{n}.csv' for n in range(12,25)]
dic_urgencias = ['/DEIS/DICCIONARIO_ATENCIONES_DE_URGENCIA_2019.xlsx',
                 '/DEIS/DICCIONARIO_ATENCIONES_DE_URGENCIA_2020.xlsx',
                 '/DEIS/DICCIONARIO_ATENCIONES_DE_URGENCIA_2021.xlsx',
                 '/DEIS/DICCIONARIO_ATENCIONES_DE_URGENCIA_2022.xlsx']
serieA = ['/DEIS/SerieA_2019.csv',
          '/DEIS/SerieA_2020.csv',
          '/DEIS/SerieA_2021.csv',
          '/DEIS/SerieA_2022.csv',
          '/DEIS/SerieA_2023.csv']

#Output
egresos_full = '/Egresos/egresos_full.parquet'
egresos_48_meses_nirse = '/Egresos/egresos_48_meses_nirse.parquet' 
egresos_48_meses = '/Egresos/egresos_48_meses.parquet' 
egresos_1461_48 = '/Egresos/egresos_1461_48.parquet'
egresos_1461_dias = '/Egresos/egresos_1461_dias.parquet' 
egresos_respiratory_causes = '/Egresos/egresos_respiratory_causes.parquet' 
columnas_original = '/columnas_egresos.csv' 
comunas_dic = '/comunas_dic.csv' 
estab_loc = '/estab_loc.csv' 
filtros = '/filtros.csv' 
trib = '/trib.csv'



#projDir = os.path.abspath(os.path.join(current_working_directory, ''))
#dataDir = current_working_directory+'/Data'

#if local_path not in sys.path:
#    sys.path.append(local_path)
if projDir not in sys.path:
    sys.path.append(projDir)
if dataDir not in sys.path:
    sys.path.append(dataDir)



#Timezones
utc = tz.gettz('UTC')
localtz = tz.gettz('America/Santiago')
today_utc = datetime.datetime.now().replace(tzinfo=utc)

today = today_utc.astimezone(localtz)
now_str = today.strftime("%Y-%m-%d %H:%M:%S")
today_str = today.strftime("%Y-%m-%d")
today_hour = today.hour


# DICCIONARIOS, LISTAS Y DATAFRAMES ÚTILES
currentWeek = today_utc.today().isocalendar().week
dropWeek = currentWeek - 1

## Listado de Causas respiratorias
# Ver como dejar esto como un archivo, no me gusta dejar en bruto
respiratorio = ["CAUSAS SISTEMA RESPIRATORIO", "TOTAL CAUSAS SISTEMA RESPIRATORIO", 'Covid-19, Virus identificado U07.1', 'Covid-19, Virus no identificado U07.2', ' - COVID-19, VIRUS IDENTIFICADO U07.1', ' - COVID-19, VIRUS NO IDENTIFICADO U07.2', 'COVID 19 Confirmado (U07.1)', 'COVID 19 Sospechoso (U07.2)', ' - COVID 19 CONFIRMADO (U07.1)', ' - COVID 19 SOSPECHOSO (U07.2)']
causasResp = ["Otra causa respiratoria (J22, J30-J39, J47, J60-J98)", "Bronquitis/bronquiolitis aguda (J20-J21)", "Influenza (J09-J11)",  "Crisis obstructiva bronquial (J40-J46)", "IRA Alta (J00-J06)", "Neumonía (J12-J18)", " - COVID-19, VIRUS IDENTIFICADO U07.1", " - COVID-19, VIRUS NO IDENTIFICADO U07.2",  "Covid-19, Virus no identificado U07.2", "Covid-19, Virus identificado U07.1", ' - COVID 19 SOSPECHOSO (U07.2)', 'COVID 19 Confirmado (U07.1)', 'COVID 19 Sospechoso (U07.2)', ' - COVID 19 CONFIRMADO (U07.1)', "TOTAL CAUSAS SISTEMA RESPIRATORIO", "CAUSAS SISTEMA RESPIRATORIO"]
causasIRAG = ["Bronquitis/bronquiolitis aguda (J20-J21)", "Influenza (J09-J11)",  "Crisis obstructiva bronquial (J40-J46)", "Neumonía (J12-J18)", " - COVID-19, VIRUS IDENTIFICADO U07.1", " - COVID-19, VIRUS NO IDENTIFICADO U07.2",  "Covid-19, Virus no identificado U07.2", "Covid-19, Virus identificado U07.1", ' - COVID 19 SOSPECHOSO (U07.2)', 'COVID 19 Confirmado (U07.1)', 'COVID 19 Sospechoso (U07.2)', ' - COVID 19 CONFIRMADO (U07.1)']
dicCausas = {"Otra causa respiratoria (J22, J30-J39, J47, J60-J98)":"Otros", "Bronquitis/bronquiolitis aguda (J20-J21)":"B/B", "Influenza (J09-J11)":"Influenza",  "Crisis obstructiva bronquial (J40-J46)":"COB", "IRA Alta (J00-J06)":"IRA Alta", "Neumonía (J12-J18)":"Neumonía", " - COVID-19, VIRUS IDENTIFICADO U07.1":"COVID", " - COVID-19, VIRUS NO IDENTIFICADO U07.2":"COVID",  "Covid-19, Virus no identificado U07.2":"COVID", "Covid-19, Virus identificado U07.1":"COVID", ' - COVID 19 SOSPECHOSO (U07.2)':"COVID", 'COVID 19 Confirmado (U07.1)':"COVID", 'COVID 19 Sospechoso (U07.2)':"COVID", ' - COVID 19 CONFIRMADO (U07.1)':"COVID", "TOTAL CAUSAS SISTEMA RESPIRATORIO":"Total", "CAUSAS SISTEMA RESPIRATORIO":"Indicaciones"}
causa_mapping = {"Bronquitis/bronquiolitis aguda (J20-J21)":"B/B", "Influenza (J09-J11)":"Influenza"}

## Diccionario de servicios
dicServicios = {"Servicio de Salud Arica y Parinacota":"SSARICA", "Servicio de Salud Tarapacá":"SSTARAP",  "Servicio de Salud Antofagasta":"SSANT", "Servicio de Salud Atacama":"SSATA","Servicio de Salud Coquimbo":"SSCOQ", 'Servicio de Salud Valparaíso San Antonio':"SSVSA", 'Servicio de Salud Viña del Mar Quillota':"SSVMQ",'Servicio de Salud Aconcagua':"SSACON",'Servicio de Salud Metropolitano Norte':"SSMN",'Servicio de Salud Metropolitano Oriente':"SSMO", 'Servicio de Salud Metropolitano Sur Oriente':"SSMSO", 'Servicio de Salud Metropolitano Sur':"SSMS", 'Servicio de Salud Metropolitano Occidente':"SSMOC", 'Servicio de Salud Metropolitano Central':"SSMC", 'Servicio de Salud Del Libertador B.O\'Higgins':"SSOH", 'Servicio de Salud Del Maule':"SSMAULE", 'Servicio de Salud Ñuble':"SSNUBLE",'Servicio de Salud Concepción':"SSCONCEP", 'Servicio de Salud Talcahuano':"SSTHNO", 'Servicio de Salud Biobío':"SSBIOBIO", 'Servicio de Salud Arauco':"SSARAUCO", 'Servicio de Salud Araucanía Norte':"SSANORTE", 'Servicio de Salud Araucanía Sur':"SSASUR", 'Servicio de Salud Los Rios':"SSLR", 'Servicio de Salud Osorno':"SSOSO", 'Servicio de Salud Del Reloncaví':"SSRELON", 'Servicio de Salud Chiloé':"SSCHL", 'Servicio de Salud Aysén':"SSAYSEN", 'Servicio de Salud Magallanes':"SSMAG"}


## Servicios Y regiones ordenados
#servicios = ['SSARICA', 'SSTARAP', 'SSANT', 'SSATA', 'SSCOQ', 'SSVSA', 'SSVMQ', 'SSACON', 'SSMN', 'SSMOC', 'SSMC', 'SSMO', 'SSMS', 'SSMSO', 'SSOH', 'SSMAULE', 'SSNUBLE', 'SSCONCEP', 'SSARAUCO', 'SSTHNO', 'SSBIOBIO', 'SSANORTE', 'SSASUR', 'SSLR', 'SSOSO', 'SSRELON', 'SSCHL', 'SSAYSEN', 'SSMAG']
#regiones = ['Arica y Parinacota', 'Tarapacá', 'Antofagasta', 'Atacama', 'Coquimbo', 'Valparaíso', 'Metropolitana', "O'Higgins", 'Maule', 'Ñuble', 'Biobío', 'Araucanía', 'Los Ríos', 'Los Lagos', 'Aysén', 'Magallanes']


#inputs nirse.cl
# datos paciente
sexo = {2:'Female',1:'Male',3:'intersex',99:'desconocido'}
prevision = {'01':'FONASA','02':'ISAPRE','03':'CAPREDENA','04':'DIPRECA','05':'SISA','96':'NINGUNA','99':'DESCONOCIDO'}
condicion_egreso = {'1':'Vivo','2':'Fallecido'}
tipo_edad = {'1':'años','2':'meses','3':'días','4':'horas'}


areasUPC = [406, 412, 415, 405, 411, 414]
areasMB = [401, 402, 403, 404, 407, 408, 409, 410, 413]
areasUPC2 = [406, 412, 415, 405, 411, 414, 310, 311, 312, 320, 323, 324]

colsDiagnostico = ['AREA_FUNC_I','AREAF_1_TRAS', 'AREAF_2_TRAS', 'AREAF_3_TRAS', 'AREAF_4_TRAS', 'AREAF_5_TRAS', 'AREAF_6_TRAS', 'AREAF_7_TRAS', 'AREAF_8_TRAS', 'AREAF_9_TRAS']
tras_date = {'AREA_FUNC_I': 'fechaIng','AREAF_1_TRAS':'fecha_tras_1', 'AREAF_2_TRAS':'fecha_tras_2', 'AREAF_3_TRAS':'fecha_tras_3', 'AREAF_4_TRAS':'fecha_tras_4', 'AREAF_5_TRAS':'fecha_tras_5'
             , 'AREAF_6_TRAS':'fecha_tras_6', 'AREAF_7_TRAS':'fecha_tras_7', 'AREAF_8_TRAS':'fecha_tras_8', 'AREAF_9_TRAS':'fecha_tras_9'}




dict_regiones = {
'Metropolitana de Santiago': 'METROPOLITANA',
'De Los Lagos': 'LOS LAGOS',
'De Valparaíso': 'VALPARAISO',
'Extranjero': 'EXTRANJERO',
'De Tarapacá': 'TARAPACA',
'Del Maule': 'MAULE',
'De Ñuble': 'NUBLE',
'Del Bíobío': 'BIOBIO',
"Del Libertador B. O'Higgins": "O'HIGGINS",
'De La Araucanía': 'ARAUCANIA',
'De Aisén del Gral. C. Ibáñez del Campo': 'AISEN',
'De Coquimbo': 'COQUIMBO',
'De Arica y Parinacota': 'ARICA Y PARINACOTA',
'De Antofagasta': 'ANTOFAGASTA',
'De Magallanes y de La Antártica Chilena': 'MAGALLANES Y ANTARTICA',
'De Los Ríos': 'LOS RIOS',
'Ignorada' : None,
'De Atacama': 'ATACAMA'
}

dict_regiones3 = {'Región Metropolitana de Santiago':'METROPOLITANA',
                    'Región de la Araucanía':'ARAUCANIA',
                    'Región del Maule':'MAULE',
                    'Región de Los Lagos':'LOS LAGOS',
                    'Antofagasta':'ANTOFAGASTA',
                    'Valparaíso':'VALPARAISO',
                    'Región del Libertador Gral. Bernardo O?Higgins':"O'HIGGINS",
                    'Región de Ñuble':'NUBLE',
                    'Región del Biobío':'BIOBIO',
                    'Arica y Parinacota':'ARICA Y PARINACOTA',
                    'Coquimbo':'COQUIMBO',
                    'Atacama':'ATACAMA',
                    'Región Aisén del Gral. Carlos Ibáñez del Campo':'AISEN',
                    'Región de Magallanes y de la Antártica Chilena':'MAGALLANES Y ANTARTICA',
                    'Región de Los Ríos':'LOS RIOS',
                    'Tarapacá':'TARAPACA'}

region_a_macrozona = {
    "ARICA Y PARINACOTA": "Macrozona Norte",
    "TARAPACA": "Macrozona Norte",
    "ANTOFAGASTA": "Macrozona Norte",
    "ATACAMA": "Macrozona Norte",  # Atacama se cubre con Tarapacá
    "COQUIMBO": "Macrozona Centro",  # Coquimbo se cubre con Valparaíso
    "VALPARAISO": "Macrozona Centro",
    "METROPOLITANA": "Macrozona METROPOLITANA",
    "O'HIGGINS": "Macrozona Centro Sur",
    #"Del Libertador B. O'Higgins":  "Macrozona Centro Sur",
    #"Del Libertador Gral. B. O'Higgins": "Macrozona Centro Sur",
    "MAULE": "Macrozona Centro Sur",  # Maule se cubre con Biobío
    "NUBLE": "Macrozona Centro Sur",
    "BIOBIO": "Macrozona Centro Sur",
    "ARAUCANIA": "Macrozona Sur",  # Araucanía se cubre con Los Ríos y Los Lagos
    "LOS RIOS": "Macrozona Sur",
    "LOS LAGOS": "Macrozona Sur",
    "AISEN": "Macrozona Austral",
    "MAGALLANES Y ANTARTICA": "Macrozona Austral"#,
    #'De Aisén del Gral. C. Ibáñez del Campo': "Macrozona Austral",
} 

region_a_macrozona2 = {
    "ARICA Y PARINACOTA": "Norte",
    "TARAPACA": "Norte",
    "ANTOFAGASTA": "Norte",
    "ATACAMA": "Norte",  # Atacama se cubre con Tarapacá
    "COQUIMBO": "Centro",  # Coquimbo se cubre con Valparaíso
    "VALPARAISO": "Centro",
    "METROPOLITANA": "Centro",
    "O'HIGGINS": "Centro",
    #"Del Libertador B. O'Higgins":  "Macrozona Centro Sur",
    #"Del Libertador Gral. B. O'Higgins": "Macrozona Centro Sur",
    "MAULE": "Centro",  # Maule se cubre con Biobío
    "NUBLE": "Centro",
    "BIOBIO": "Centro",
    "ARAUCANIA": "Sur",  # Araucanía se cubre con Los Ríos y Los Lagos
    "LOS RIOS": "Sur",
    "LOS LAGOS": "Sur",
    "AISEN": "Sur",
    "MAGALLANES Y ANTARTICA": "Sur"#,
    #'De Aisén del Gral. C. Ibáñez del Campo': "Macrozona Austral",
} 


# Códigos agrupados por RSV
diagnosticosVRS=['J121', 'J205', 'J210','J219', 'B974' ]

rsv_0 = [
    "J121",  # Pneumonia Due to Respiratory Syncytial Virus
    "J205",  # Acute Bronchitis Due to Respiratory Syncytial Virus
    "J210"   # Acute Bronchiolitis Due to Respiratory Syncytial Virus
]

rsv_1 = rsv_0 + [
    "J219",  # Acute Bronchiolitis, Unspecified
    "B974"   # Respiratory Syncytial Virus as Cause of Diseases Classified in Other Chapters (*)
]

rsv_2 = rsv_1 + [
    "J128",  # Pneumonia Due to Other Viruses
    "J129",  # Viral Pneumonia, Unspecified (Virus)
    "J180",  # Bronchopneumonia, Unspecified
    "J188",  # Other Pneumonias, of Unspecified Microorganism
    "J189",  # Pneumonia, Unspecified
    "J208",  # Acute Bronchitis Due to Other Specified Microorganisms
    "J209",  # Acute Bronchitis, Unspecified
    "J218",  # Acute Bronchiolitis Due to Other Specified Microorganisms
    "J22X"   # Unspecified Acute Infection of Lower Respiratory Tract
]

filter_estab = [ 120103., 124130., 105107., 103102., 103203., 120105.,
       107206., 103219., 108204., 107108., 114204., 128112., 126204.,
       105106., 125103., 115206., 115221., 112211., 107223., 106205.,
       102201., 121110., 105208., 111101., 112102., 107224., 124260.,
       122202., 124210., 112212., 109200., 200486., 101213., 112243.,
       201288., 200234., 104101., 117104., 201319., 126704.]

tipo_edad_dic = {
    1: 'años',
    2: 'meses',
    3: 'días',
    4: 'horas'
} 
tipo_edad_dic_2 = {
    'años': 12,
    'meses': 1,
    'días': 1/30,
    'horas': (1/24)/30
} 

etnia_dic = {
    1: 'MAPUCHE',
    2: 'AYMARA',
    3: 'RAPA NUI (PASCUENSE)',
    4: 'LICAN ANTAI (ATACAMEÑO)',
    5: 'QUECHUA',
    6: 'COLLA',
    7: 'DIAGUITA',
    8: 'KAWÉSQAR',
    9: 'YAGÁN (YÁMANA)',
    10: 'OTRO',
    11: 'NO SABE / NO RESPONDE',
    96: 'NINGUNO'
}

trib_servicios_dic ={'Arica y Parinacota': 'SSARICA',
                        'Arica': 'SSARICA',
                        'Tarapacá': 'SSTARAP',
                        'Antofagasta': 'SSANT',
                        'Atacama': 'SSATA',
                        'Coquimbo': 'SSCOQ',
                        'Valparaíso': 'SSVSA',
                        'Metropolitana Santiago':  '[SSMN, SSMO, SSMSO, SSMS, SSMOC, SSMC]',
                        'del Libertador B. O Higgins': 'SSOH',
                        'del Maule': 'SSMAULE',
                        'Ñuble': 'SSNUBLE',
                        'del Biobío': 'SSBIOBIO',
                        'La Araucanía': '[SSANORTE, SSASUR]',
                        'Los Rios': 'SSLR',
                        'Los Lagos': '[SSOSO, SSRELON, SSCHL]',
                        'Magallanes y la Antártica Chilena': 'SSMAG',
                        'Valparaíso San Antonio': 'SSVSA',
                        'Viña del Mar Quillota': 'SSVMQ',
                        'Aconcagua': 'SSACON',
                        'Metropolitano Norte': 'SSMN',
                        'Metropolitano Occidente': 'SSMOC',
                        'Metropolitano Central': 'SSMC',
                        'Metropolitano Oriente': 'SSMO',
                        'Metropolitano Sur': 'SSMS',
                        'Metropolitano Sur Oriente': 'SSMSO',
                        'Del Libertador B.O Higgins': 'SSOH',
                        'Del Maule': 'SSMAULE',
                        'Concepción': 'SSCONCEP',
                        'Arauco': 'SSARAUCO',
                        'Talcahuano': 'SSTHNO',
                        'Biobío': 'SSBIOBIO',
                        'Araucanía Norte': 'SSANORTE',
                        'Araucanía Sur': 'SSASUR',
                        'Osorno': 'SSOSO',
                        'Del Reloncaví': 'SSRELON',
                        'Chiloé': 'SSCHL',
                        'Aisén': 'SSAYSEN',
                        'Magallanes': 'SSMAG'}
## Diccionario de servios a regiones
trib_servicios_reg_dic = {'SSARICA': 'Arica y Parinacota',
                          'SSARICA': 'Arica y Parinacota',
                          'SSTARAP': 'Tarapacá',
                          'SSANT': 'Antofagasta',
                          'SSATA': 'Atacama',
                          'SSCOQ': 'Coquimbo',
                          'SSVSA': 'Valparaíso',
                          '[SSMN, SSMO, SSMSO, SSMS, SSMOC, SSMC]': 'Metropolitana',
                          'SSOH': 'OHiggins',
                          'SSMAULE': 'Maule',
                          'SSNUBLE': 'Ñuble',
                          'SSBIOBIO': 'Biobío',
                          '[SSANORTE, SSASUR]': 'Araucanía',
                          'SSLR': 'Los Ríos',
                          '[SSOSO, SSRELON, SSCHL]': 'Los Lagos',
                          'SSMAG': 'Magallanes',
                          'SSVSA': 'Valparaíso',
                          'SSVMQ': 'Valparaíso',
                          'SSACON': 'Valparaíso',
                          'SSMN': 'Metropolitana',
                          'SSMOC': 'Metropolitana',
                          'SSMC': 'Metropolitana',
                          'SSMO': 'Metropolitana',
                          'SSMS': 'Metropolitana',
                          'SSMSO': 'Metropolitana',
                          'SSOH': 'OHiggins',
                          'SSMAULE': 'Maule',
                          'SSCONCEP': 'Biobío',
                          'SSARAUCO': 'Biobío',
                          'SSTHNO': 'Biobío',
                          'SSBIOBIO': 'Biobío',
                          'SSANORTE': 'Araucanía',
                          'SSASUR': 'Araucanía',
                          'SSOSO': 'Los Lagos',
                          'SSRELON': 'Los Lagos',
                          'SSCHL': 'Los Lagos',
                          'SSAYSEN': 'Aysén',
                          'SSMAG': 'Magallanes'}

region_a_macrozona_V2 = {
    "Arica y Parinacota": "Macrozona Norte",
    "Tarapacá": "Macrozona Norte",
    "Antofagasta": "Macrozona Norte",
    "Atacama": "Macrozona Norte",  # Atacama se cubre con Tarapacá
    "Coquimbo": "Macrozona Centro",  # Coquimbo se cubre con Valparaíso
    "Valparaíso": "Macrozona Centro",
    "Metropolitana": "Macrozona METROPOLITANA",
    "OHiggins": "Macrozona Centro Sur",
    "Maule": "Macrozona Centro Sur",  # Maule se cubre con Biobío
    "Ñuble": "Macrozona Centro Sur",
    "Biobío": "Macrozona Centro Sur",
    "Araucanía": "Macrozona Sur",  # Araucanía se cubre con Los Ríos y Los Lagos
    "Los Ríos": "Macrozona Sur",
    "Los Lagos": "Macrozona Sur",
    "Aysén": "Macrozona Austral",
    "Magallanes": "Macrozona Austral"
} 

region_a_macrozona2_V2 = {
    "Arica y Parinacota": "Norte",
    "Tarapacá": "Norte",
    "Antofagasta": "Norte",
    "Atacama": "Norte",  # Atacama se cubre con Tarapacá
    "Coquimbo": "Centro",  # Coquimbo se cubre con Valparaíso
    "Valparaíso": "Centro",
    "Metropolitana": "Centro",
    "OHiggins": "Centro",
    "Maule": "Centro",  # Maule se cubre con Biobío
    "Ñuble": "Centro",
    "Biobío": "Centro",
    "Araucanía": "Sur",  # Araucanía se cubre con Los Ríos y Los Lagos
    "Los Ríos": "Sur",
    "Los Lagos": "Sur",
    "Aysén": "Sur",
    "Magallanes": "Sur"#
} 

dict_regiones2 = {'Región Metropolitana de Santiago':'Metropolitana',
                  "Región Del Libertador Gral. B. O'Higgins":'OHiggins',
                  'Región De Ñuble':'Ñuble',
                  'Región Del Bíobío':'Biobío',
                  'Región De Los Ríos':'Los Ríos',
                  'Región De Los Lagos':'Los Lagos',
                  'Región De Aysén del General Carlos Ibañez del Campo': 'Aysén',
                  'Región De La Araucanía':'Araucanía',
                  'Región De Arica Parinacota':'Arica y Parinacota',
                  'Región De Antofagasta':'Antofagasta',
                  'Región De Atacama':'Atacama',
                  'Región De Coquimbo':'Coquimbo',
                  'Región De Valparaíso':'Valparaíso',
                  'Región De Magallanes y de la Antártica Chilena':'Magallanes',
                  'Región De Tarapacá':'Tarapacá',
                  'Región Del Maule':'Maule'}

cat_ICD10 = ['Acute Bronchitis/Bronchiolitis',
                'Pneumonia',
                'Influenza',
                'Upper respiratory tract infection',
                'Obtructive Pulmonary Disease',
                'COVID-19',
                'Other respiratory conditions']

urgencia_causas = [3,              # Acute Bronchitis/Bronchiolitis
                    5,              # Pneumonia
                    4,              # Influenza
                    10,             # Upper respiratory tract infection
                    11,             # Obtructive Pulmonary Disease
                    30,31,32,33,    # COVID-19
                    6]              # Other respiratory conditions

urgencias_causas_dic = {3: 'Acute Bronchitis/Bronchiolitis',
                        5: 'Pneumonia',
                        4: 'Influenza',
                        10: 'Upper respiratory tract infection',
                        11: 'Obtructive Pulmonary Disease',
                        30: 'COVID-19',
                        31: 'COVID-19',
                        32: 'COVID-19',
                        33: 'COVID-19',
                        6: 'Other respiratory conditions'}

REM_serie_A_consultas = ['03020101',   #IRA ALTA
                        '03020201',  #SÍNDROME BRONQUIAL OBSTRUCTIVO
                        '03020301',	 #NEUMONÍA
                        '03020402',	 #ASMA
                        '03020403',	 #ENFERMEDAD PULMONAR OBSTRUCTIVA CRÓNICA
                        '03020401',  #OTRAS RESPIRATORIAS
                        '03040210',	 #OBSTÉTRICA
                        '03040220',	 #GINECOLÓGICA
                        '04040100',	 #GINECOLÓGICA  POR INFERTILIDAD
                        '04025010',	 #INFECCIÓN TRANSMISIÓN SEXUAL
                        '04025020',	 #VIH-SIDA
                        '04025025',	 #SALUD MENTAL
                        '04040427',	 #CARDIOVASCULAR
                        '03020501']	 #OTRAS MORBILIDADES

REM_serie_A_consultas_cat_ICD10 = ['03020101',
                                    '03020201',
                                    '03020301',
                                    '03020402',
                                    '03020403',
                                    '03020401']
dic_REM_cat_ICD10 = {'03020101': 'Upper respiratory tract infection',
                        '03020201': 'Acute Bronchitis/Bronchiolitis',
                        '03020301':	'Pneumonia',
                        '03020402':	'Other respiratory conditions',
                        '03020403':	'Acute Bronchitis/Bronchiolitis',
                        '03020401':	'Other respiratory conditions'}

# REM 2019-2022
REM_serie_A_nacimientos_peso = ['01060100',	 #NACIDOS VIVOS    
                                '01060101']  #NACIDOS FALLECIDOS
REM_serie_A_nacimientos_malformacion = ['24200110',	 #NACIDOS VIVOS    
                                        '24200120']  #NACIDOS FALLECIDOS
REM_serie_A_nacimientos_atendidos = ['01030100',	 #NORMAL/VAGINAL
                                    '01030200',	     #DISTÓCICO VAGINAL
                                    '01030300', 	 #CESÁREA ELECTIVA
                                    '24090700',	     #CESÁREA URGENCIA
                                    '01030400',	     #ABORTOS
                                    '24090800',	     #PARTO NORMAL VERTICAL
                                    '24200121',	     #ENTREGA DE PLACENTA A SOLICITUD DE LA MUJER
                                    '24090900',	     #PARTO FUERA ESTABLECIMIENTO DE SALUD
                                    '24100100']	     #EMBARAZO NO CONTROLADO
dic_REM_serie_A = {'03020101': 	'IRA ALTA',
                    '03020201': 'SÍNDROME BRONQUIAL OBSTRUCTIVO',
                    '03020301':	'NEUMONÍA',
                    '03020402':	'ASMA',
                    '03020403':	'ENFERMEDAD PULMONAR OBSTRUCTIVA CRÓNICA',
                    '03020401':	'OTRAS RESPIRATORIAS',
                    '03040210':	'OBSTÉTRICA',
                    '03040220':	'GINECOLÓGICA',
                    '04040100':	'GINECOLÓGICA  POR INFERTILIDAD',
                    '04025010':	'INFECCIÓN TRANSMISIÓN SEXUAL',
                    '04025020':	'VIH-SIDA',
                    '04025025':	'SALUD MENTAL',
                    '04040427':	'CARDIOVASCULAR',
                    '03020501':	'OTRAS MORBILIDADES',
                    '01060100':	'NACIDOS VIVOS',       
                    '01060101':	'NACIDOS FALLECIDOS',
                    '24200110': 'NACIDOS VIVOS',
                    '24200120': 'NACIDOS FALLECIDOS',
                    '01030100':	'NORMAL/VAGINAL',      
                    '01030200':	'DISTÓCICO VAGINAL',
                    '01030300':	'CESÁREA ELECTIVA',
                    '24090700':	'CESÁREA URGENCIA',
                    '01030400': 'ABORTOS',
                    '24090800': 'PARTO NORMAL VERTICAL',
                    '24200121':	'ENTREGA DE PLACENTA A SOLICITUD DE LA MUJER',
                    '24090900':	'PARTO FUERA ESTABLECIMIENTO DE SALUD',
                    '24100100':	'EMBARAZO NO CONTROLADO'}

# REM 2023
REM_serie_A_nacimientos_peso_23 = ['01060100']	 #NACIDOS VIVOS    
		
REM_serie_A_nacimientos_atendidos_23 = ['01030100',	     #NORMAL/VAGINAL
                                        '01030300', 	 #CESÁREA ELECTIVA
                                        '24090700',	     #CESÁREA URGENCIA
                                        '24200121',	     #ENTREGA DE PLACENTA A SOLICITUD DE LA MUJER
                                        '24100100',      #EMBARAZO NO CONTROLADO
                                        '29101714',      #INSTRUMENTAL
                                        '29101715',      #PARTO PREHOSPITALARIO (en establecimientos salud o ambulancias)	
                                        '29101716',      #PARTOS FUERA DE LA RED DE SALUD	
                                        '29101717',      #PARTO EN DOMICILIO CON ATENCION PROFESIONAL
                                        '29101718']      #PARTO EN DOMICILIO SIN ATENCION PROFESIONAL
                                           
dic_REM_serie_A_2023 = {'01060100':	'NACIDOS VIVOS',       
                        '01030100':	'NORMAL/VAGINAL',      
                        '01030300':	'CESÁREA ELECTIVA',
                        '24090700':	'CESÁREA URGENCIA',
                        '24200121':	'ENTREGA DE PLACENTA A SOLICITUD DE LA MUJER',
                        '24100100':	'EMBARAZO NO CONTROLADO',
                        '29101714': 'INSTRUMENTAL',
                        '29101715': 'PARTO PREHOSPITALARIO',
                        '29101716': 'PARTOS FUERA DE LA RED DE SALUD',
                        '29101717': 'PARTO EN DOMICILIO CON ATENCION PROFESIONAL',
                        '29101718': 'PARTO EN DOMICILIO SIN ATENCION PROFESIONAL'}



min_fecha = pd.to_datetime({'year': [2010], 'month': [1], 'day': [1]}, format='%Y-%m-%d')[0]
max_fecha = pd.to_datetime({'year': [2029], 'month': [12], 'day': [31]}, format='%Y-%m-%d')[0]


p_basic_bed = 414
p_icu_bed = 1083
emergency = 192
nirsevimab = 225
dose_administration = 5 #only catch-up
palivizumab_19 = 10683035
palivizumab_22 = 11321037
palivizumab_23 = 12568198


cardiopatias = [
    'Q200', 'Q201', 'Q202', 'Q203', 'Q204', 'Q205', 'Q206', 'Q208', 'Q209',
    'Q210', 'Q211', 'Q212', 'Q213', 'Q214', 'Q218', 'Q219',
    'Q220', 'Q221', 'Q222', 'Q223', 'Q224', 'Q225', 'Q226', 'Q228', 'Q229',
    'Q230', 'Q231', 'Q232', 'Q233', 'Q234', 'Q238', 'Q239',
    'Q240', 'Q242', 'Q243', 'Q244', 'Q245', 'Q248', 'Q249',
    'Q250', 'Q251', 'Q252', 'Q253', 'Q254', 'Q255', 'Q256', 'Q257', 'Q258', 'Q259',
    'Q260', 'Q262', 'Q263', 'Q264'
]

cardiopatias2 = [
    'Q200', 'Q201', 'Q202', 'Q203', 'Q204', 'Q205', 'Q206', 'Q208', 'Q209',
    'Q210', 'Q212', 'Q213', 'Q214', 'Q218', 'Q219',
    'Q220', 'Q221', 'Q222', 'Q223', 'Q224', 'Q225', 'Q226', 'Q228', 'Q229',
    'Q230', 'Q231', 'Q232', 'Q234', 'Q238', 'Q239',
    'Q240', 'Q242', 'Q243', 'Q244', 'Q245', 'Q248', 'Q249',
    'Q251', 'Q252', 'Q253', 'Q254', 'Q255', 'Q256', 'Q257', 'Q258', 'Q259',
    'Q260', 'Q262', 'Q263', 'Q264'
]

displacia_pulmonar = ['P271']



