import pandas as pd
import os
from pathlib import Path
import numpy as np
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
from difflib import SequenceMatcher

def son_similares(cadena1, cadena2, umbral=0.6):
    similitud = SequenceMatcher(None, cadena1, cadena2).ratio()
    return similitud >= umbral

path_actual = Path.cwd()
path_data = path_actual/'Nirse_cl' / 'Data'
print("PATH Data:",path_data)
new = pd.read_csv(path_data / '2024-08-27_Egresos.csv',encoding= 'latin1', sep = '|')
egresos = pd.read_csv(path_data / "egresos.csv",encoding = "latin1",sep="|")

#egresos = egresos[egresos['ANO_EGR']!=2024]
egresos = pd.concat([egresos,new])

trib_publi = pd.read_excel(path_data/"tributacion.xlsx")
trib_priv = pd.read_excel(path_data/"EstadoCargaIEEH_SEREMI_19082024.xlsx", skiprows=2)
trib = pd.concat([trib_publi,trib_priv])
trib['CodigoEstablecimiento'] = trib['Codigo Establecimiento']

acusados = trib[trib['Jun'].isna()].CodigoEstablecimiento
acusados = acusados.dropna().astype(int)

diagnosticos_upc = [406, 412, 415, 405, 411, 414]
cols_diagnostico = ['AREA_FUNC_I','AREAF_1_TRAS', 'AREAF_2_TRAS', 'AREAF_3_TRAS', 'AREAF_4_TRAS', 'AREAF_5_TRAS', 'AREAF_6_TRAS', 'AREAF_7_TRAS', 'AREAF_8_TRAS', 'AREAF_9_TRAS']
tras_date = {'AREA_FUNC_I': 'fechaIng','AREAF_1_TRAS':'fecha_tras_1', 'AREAF_2_TRAS':'fecha_tras_2', 'AREAF_3_TRAS':'fecha_tras_3', 'AREAF_4_TRAS':'fecha_tras_4', 'AREAF_5_TRAS':'fecha_tras_5'
             , 'AREAF_6_TRAS':'fecha_tras_6', 'AREAF_7_TRAS':'fecha_tras_7', 'AREAF_8_TRAS':'fecha_tras_8', 'AREAF_9_TRAS':'fecha_tras_9'}

def obtener_fecha_primer_upc(row):
    for col in cols_diagnostico:
        if row[col] == 1:
            fecha_col = tras_date[col]
            return row[fecha_col]
    return None

egresos['FECHA_NAC'] = pd.to_datetime({'year': egresos['A_NAC'], 'month': egresos['M_NAC'], 'day': egresos['D_NAC']}, format='%Y-%m-%d')
egresos['fechaIng'] = pd.to_datetime({'year': egresos['ANO_ING'], 'month': egresos['MES_ING'], 'day': egresos['DIA_ING']}, format='%Y-%m-%d')

#e24 = e24[e24.fechaIng <= "2024-04-21"]
elementos = ['J121', 'J205', 'J210','J219', 'B974' ]
df_filtrado = egresos[egresos[['DIAG1'#,'DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9','DIAG10','DIAG11'
                               ]].isin(elementos).any(axis=1)]

df_filtrado['VRS1'] = 1
df_filtrado=df_filtrado.rename(columns={'RUT':'RUN'})

comunas = pd.read_excel(path_data/"comunas.xlsx")
comunas = comunas.drop_duplicates()
comunas = comunas.rename(columns = {'C_COM': 'COMUNA','NOM_REG':'NOMBRE_REGION'})
df_filtrado = df_filtrado.merge(comunas,how='left',on ='COMUNA')
df_filtrado['NOMBRE_REGION'] = df_filtrado['NOMBRE_REGION'].str.replace('Región', '').str.strip()
df_filtrado['NOMBRE_REGION'] = df_filtrado['NOMBRE_REGION'].fillna('Ignorada')
df_filtrado = df_filtrado[['RUN', 'FECHA_NAC', 'fechaIng', 'VRS1','SEXO','NOMBRE_REGION','AREA_FUNC_I', 'ANO_ING','ESTAB',
                           'DIA_1_TRAS', 'MES_1_TRAS', 'ANO_1_TRAS', 'AREAF_1_TRAS',
                            'DIA_2_TRAS', 'MES_2_TRAS', 'ANO_2_TRAS', 'AREAF_2_TRAS', 'DIA_3_TRAS',
                            'MES_3_TRAS', 'ANO_3_TRAS', 'AREAF_3_TRAS', 'DIA_4_TRAS', 'MES_4_TRAS',
                            'ANO_4_TRAS', 'AREAF_4_TRAS', 'DIA_5_TRAS', 'MES_5_TRAS', 'ANO_5_TRAS',
                            'AREAF_5_TRAS', 'DIA_6_TRAS', 'MES_6_TRAS', 'ANO_6_TRAS',
                            'AREAF_6_TRAS', 'DIA_7_TRAS', 'MES_7_TRAS', 'ANO_7_TRAS',
                            'AREAF_7_TRAS', 'DIA_8_TRAS', 'MES_8_TRAS', 'ANO_8_TRAS',
                            'AREAF_8_TRAS', 'DIA_9_TRAS', 'MES_9_TRAS', 'ANO_9_TRAS',
                            'AREAF_9_TRAS']]

df_filtrado=df_filtrado.assign(age=lambda x: ((x['fechaIng'] - x['FECHA_NAC']).dt.days / 30))
df_filtrado=df_filtrado.query("age <= 48")  ###############################################################################################################################

e24 = df_filtrado[df_filtrado.ANO_ING ==2024]
dataprev = df_filtrado[df_filtrado.ANO_ING !=2024]

e24.to_csv(path_data/"egresos2024.csv")
dataprev.to_csv(path_data/"dataprev.csv")
df_filtrado.to_csv(path_data/"full_data.csv")


data = (
    pd.read_csv(path_data/"full_data.csv")[['RUN', 'FECHA_NAC', 'fechaIng','age']]
    .assign(FECHA_NAC=lambda data: pd.to_datetime(data["FECHA_NAC"], format="%Y-%m-%d"))
    .assign(fechaIng=lambda data: pd.to_datetime(data["fechaIng"], format="%Y-%m-%d"))
    .assign(epiweek=lambda x: x['fechaIng'].dt.isocalendar().week)
    .assign(year=lambda x: x['fechaIng'].dt.isocalendar().year)
)

data = data.drop_duplicates(subset=['RUN','fechaIng'], keep='first')

catchup24 =( 
    data.query("'2023-10-01' <= FECHA_NAC <= '2024-03-31' ")
    .assign(season = lambda x: 'pre_season').assign(elegibilidad = lambda x: 2024)
)
inseason24 = (
    data.query("'2024-04-01' <= FECHA_NAC <= '2024-09-30'")
    .assign(season = lambda x: 'in_season').assign(elegibilidad = lambda x: 2024)
)

catchup23 =( 
    data.query("'2022-10-01' <= FECHA_NAC <= '2023-03-31' ")
    .assign(season = lambda x: 'pre_season').assign(elegibilidad = lambda x: 2023)
)
inseason23 = (
    data.query("'2023-04-01' <= FECHA_NAC <= '2023-09-30'")
    .assign(season = lambda x: 'in_season').assign(elegibilidad = lambda x: 2023)
)
catchup22 =( 
    data.query("'2021-10-01' <= FECHA_NAC <= '2022-03-31'")
    .assign(season = lambda x: 'pre_season').assign(elegibilidad = lambda x: 2022)
)
inseason22 = (
    data.query("'2022-04-01' <= FECHA_NAC <= '2022-09-30'")
    .assign(season = lambda x: 'in_season').assign(elegibilidad = lambda x: 2022)
)
catchup19=( 
    data.query("'2018-10-01' <= FECHA_NAC <= '2019-03-31' ")
    .assign(season = lambda x: 'pre_season').assign(elegibilidad = lambda x: 2019)
)
inseason19 = (
    data.query("'2019-04-01' <= FECHA_NAC <= '2019-09-30'")
    .assign(season = lambda x: 'in_season').assign(elegibilidad = lambda x: 2019)
)

#merge seasons
seasons = pd.concat([catchup24,inseason24,catchup23,inseason23,catchup22,inseason22,catchup19,inseason19])[['RUN','fechaIng','season','elegibilidad']]
df = data.merge(seasons, how='left', on=['RUN','fechaIng'])
df['elegibilidad'] = df['elegibilidad'].fillna(0)
df['season'] = df['season'].fillna('nonelegible')

#merge egresos
egresosprev = pd.read_csv(path_data/'dataprev.csv')
egresos24 = pd.read_csv(path_data/"egresos2024.csv")

egresosprev['fechaIng'] = pd.to_datetime(egresosprev['fechaIng'], format='%Y-%m-%d')
egresos24['fechaIng'] = pd.to_datetime(egresos24['fechaIng'], format='%Y-%m-%d')

for i in range(1, 10):
    year_col = f'ANO_{i}_TRAS'
    month_col = f'MES_{i}_TRAS'
    day_col = f'DIA_{i}_TRAS'
    date_col = f'fecha_tras_{i}'

    # Combinar las columnas de año, mes y día en una columna de fecha
    egresosprev[date_col] = pd.to_datetime({'year': egresosprev[year_col], 'month': egresosprev[month_col], 'day': egresosprev[day_col]}, format='%Y-%m-%d')
    egresos24[date_col] = pd.to_datetime({'year': egresos24[year_col], 'month': egresos24[month_col], 'day': egresos24[day_col]}, format='%Y-%m-%d')

for col in cols_diagnostico:
    egresosprev[col] = egresosprev[col].apply(lambda x: 1 if x in diagnosticos_upc else 0)
    egresos24[col] = egresos24[col].apply(lambda x: 1 if x in diagnosticos_upc else 0)

egresosprev['cama'] = np.where(egresosprev[cols_diagnostico].eq(1).any(axis=1), 'UPC', "")
egresosprev['critico'] = np.where(egresosprev[cols_diagnostico].eq(1).any(axis=1), 1, 0)
egresosprev['ingreso_UPC'] = egresosprev.apply(obtener_fecha_primer_upc, axis=1)
egresos24['critico'] = np.where(egresos24[cols_diagnostico].eq(1).any(axis=1), 1, 0)
egresos24['cama'] = np.where(egresos24[cols_diagnostico].eq(1).any(axis=1), 'UPC', "")
egresos24['ingreso_UPC'] = egresos24.apply(obtener_fecha_primer_upc, axis=1)


egresos = pd.concat([egresosprev, egresos24], ignore_index=True)

# Aplicar la transformación
egresos = egresos.drop_duplicates(subset=['RUN','fechaIng'], keep='first')

egresos = egresos[['RUN','NOMBRE_REGION','SEXO','cama','ingreso_UPC','critico','fechaIng','ESTAB']]
df = df.merge(egresos, how = 'left', on = ['RUN','fechaIng'])

regiones = {
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

sex = {'2':'Female','1':'Male','9':'intersex'}
df['SEXO'] = df['SEXO'].astype(int)
df['SEXO'] = df['SEXO'].astype(str)

#df.loc[:, 'NOMBRE_REGION'] = df['NOMBRE_REGION'].map(regiones)

def mapear_region(region, regiones_dict, umbral=0.6):
    for key, value in regiones_dict.items():
        if son_similares(region, key, umbral):
            return value
    return None 

df['NOMBRE_REGION'] = df['NOMBRE_REGION'].apply(lambda x: mapear_region(x, regiones))

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

df['Macrozona'] = df['NOMBRE_REGION'].map(region_a_macrozona)

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

df['Macrozona2'] = df['NOMBRE_REGION'].map(region_a_macrozona2)

df.loc[:, 'SEXO'] = df['SEXO'].map(sex)
df['age'] = df['age'].astype(int)

df = df.rename(columns ={'NOMBRE_REGION':'region','SEXO':'sex'})
df=df.dropna(subset=['region','sex','elegibilidad'])
df["epiweekupc"] = df['ingreso_UPC'].dt.isocalendar().week.fillna(-1).astype(int)
df = df.drop_duplicates()
df['VRS1'] = 1

df_acusados = df[df['ESTAB'].isin(acusados)] 

trib['trib_suf'] = trib['Jun'] <= 0.8 * trib['May']
df_trib = trib[~trib['CodigoEstablecimiento'].isna()]
acusados_2 = df_trib[df_trib.trib_suf].CodigoEstablecimiento

df_all = df.copy()

df = df[~df['ESTAB'].isin(acusados)]
df = df[~df['ESTAB'].isin(acusados_2)]

cortes = {17:'2024-04-21',18:'2024-04-28',19:'2024-05-05',20:'2024-05-12',21:'2024-05-19',22:'2024-05-26',23:'2024-06-02',24:'2024-06-09',25:'2024-06-16',26:'2024-06-23',27:'2024-06-30',
          28:'2024-07-07',29:'2024-07-14',30:'2024-07-21'}

df = df[df.fechaIng <= cortes[26]] #FILTRADO FECHAS PEDIDO PROFES

df.to_csv(path_data/"data.csv")

df[df.year!=2018].groupby(['year','epiweek','VRS1','elegibilidad','critico'],as_index=False)[['RUN']].count().rename(columns={'RUN':'Count'}).to_csv(path_data/'data_download.csv')


mshare_acusados = []
u=0
inter=0
max_estab=0
d_2023_25=df_all[(df_all.year==2023) #& (df_all.epiweek<=25)
                 ]
ntot = d_2023_25.shape[0]
for i in acusados.unique():
    nrow = d_2023_25[d_2023_25['ESTAB'] == i].shape[0]
    if nrow==0:
        pass
    elif u <= nrow/ntot:
        u=nrow/ntot
        max_estab = i
        mshare_acusados.append(nrow/ntot)
    else:
        mshare_acusados.append(nrow/ntot)

for i in acusados_2.unique():
    if i in acusados.unique():
        print("wow")
        inter=1
    else:
        nrow = d_2023_25[d_2023_25['ESTAB'] == i].shape[0]
        if nrow==0:
            pass
        elif u <= nrow/ntot:
            u=nrow/ntot
            max_estab = i
            mshare_acusados.append(nrow/ntot)
        else:
            mshare_acusados.append(nrow/ntot)

print("Fueron eliminados el", round(sum(mshare_acusados)*100,2), "% de los datos")
if inter==0:
    print("Fueron eliminados", acusados.nunique()+acusados_2.nunique(), "hospitales")