import pandas as pd
import datetime as dt
import numpy as np
import argparse
import sys
import warnings
import os
import datetime 
from dateutil import tz
import pandas as pd
from dateutil.relativedelta import relativedelta

# Suppress FutureWarnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

#Paths for python 
current_working_directory = os.getcwd()
projDir = os.path.abspath(os.path.join(current_working_directory, '../..', ''))
dataDir = os.path.abspath(os.path.join(current_working_directory, '..', 'Data'))
dataDirI = dataDir+'/Input'
dataDirO = dataDir+'/Output'

############################################# Info ##############################################
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

############################################# Paths ##############################################
path_info_paises = '/Ficha Activo Información Egresos_Hospitalarios 2.0.xlsx'
path_comunas = 'comunas_dic.csv'

########################################## Diccionarios ##########################################
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
sexo = {2:'Female',1:'Male',3:'intersex',99:'desconocido'}
tipo_edad_dic = {1: 'años',
                 2: 'meses',
                 3: 'días',
                 4: 'horas'} 

###################################################################################################
############################################# EGRESOS #############################################
###################################################################################################
 
###################################### Funciones auxiliares ####################################### 

# Función que lee el archivo de egresos, modifica información básica y lo transforma a parquet
def pre_egresos(path):
    egresos_all = pd.read_csv(f'{dataDirI}/{path}',encoding= 'latin1', sep = '|', low_memory=False)

    df = (egresos_all
          .rename(columns={'RUT':'RUN'})
          .assign(fechaNac =  lambda df: pd.to_datetime({'year': df['A_NAC'], 'month': df['M_NAC'], 'day': df['D_NAC']}, format='%Y-%m-%d'),
                  fechaIng = lambda df: pd.to_datetime({'year': df['ANO_ING'], 'month': df['MES_ING'], 'day': df['DIA_ING']}, format='%Y-%m-%d'))
          .assign(dias_vida =lambda x: ((x['fechaIng'] - x['fechaNac']).dt.days))
          .assign( 
                fechaEgr = lambda df: pd.to_datetime({'year': df['ANO_EGR'], 'month': df['MES_EGR'], 'day': df['DIA_EGR']}, format='%Y-%m-%d'),
                duracion = lambda df: (df['fechaEgr']-df['fechaIng']).dt.days.apply(lambda x: max(1, x)),
                epiweek=lambda x: x['fechaIng'].dt.isocalendar().week,
                year=lambda x: x['fechaIng'].dt.isocalendar().year,
                first_letter =  lambda x: x.DIAG1.str[0],
                diag_respiratory_causes=lambda x: x['DIAG1'].str.startswith('J') | x['DIAG1'].eq('B974')))
    df.to_parquet(f'{dataDirI}/egresos_sin_filtro.parquet', engine='pyarrow', compression='snappy')
    df = pd.read_parquet(f'{dataDirO}/egresos_sin_filtro.parquet', engine='pyarrow')
    return df

# Función que realiza los filtros a la base de egresos
def filtros_egresos(df):
    df = df[(df['year'] >= 2018) & 
        (~df['ESTAB'].isna()) & 
        (df['dias_vida'] >= 0) & 
        (df['DIAS_ESTAD'] >= 0) &
        (~df['fechaIng'].isna())]
    return df

# Función que agrega la información de los establecimientos
def estab(df):
    trib = pd.read_csv(f'{dataDirO}/trib.csv')
    estab_loc = (pd.read_csv(f'{dataDirO}/estab_loc.csv')[['ESTAB','COMUNA_ESTAB','comuna_estab_str','Region']].rename(columns={'Region':'Region_ESTAB'}))

    df = (df
          .merge(trib[['ESTAB',
                       'NombreEstablecimiento',
                       'origen_estab',
                       'NombreServicioSalud',
                       'NombreSeremiSalud',
                       'Region_SS',
                       'Macrozona1_SS',
                       'Macrozona2_SS',
                       'trib_rech']],how='left',on='ESTAB')
          .fillna({'trib_rech': True}))
    
    df = (df.merge(estab_loc,how='left',on='ESTAB'))
    return df

# Indicador de si el paciente nació en la temporada de invierno
def inseason(month_born,year_born,year_in):
    return (((year_in-year_born)==0)*(month_born>=4))*(month_born<=9)*1

# Función que transforma la información importante contenida en el dataset
def categorical_info(df,path_paises=path_info_paises,path_comunas=path_comunas,sexo=sexo,tipo_edad_dic=tipo_edad_dic,etnia_dic=etnia_dic):
    # Columnas a convertir
    serv_cl_col = list([f'SERC_{x}_TRAS' for x in range(1,10)])+['SER_CLIN_I','SERC_EGR']
    areaf_col = list([f'AREAF_{x}_TRAS' for x in range(1,10)])+['AREAF_1_TRAS','AREAF_EGR']
    columnas_a_convertir = ['ESTAB','ServicioSalud','Seremi','SEXO','TIPO_EDAD','ETNIA','P_ORIGEN','COMUNA','PREVI','COND_EGR']
    columnas_a_convertir2 = ['sexo_str','tipo_edad_str','etnia_str','pais_origen','nacional','COMUNA']

    # Cargamos los diccionarios
    paises_dic = pd.read_excel(f'{dataDirI}/{path_paises}', sheet_name='Anexo 3', skiprows=6)
    paises_dic = dict(zip(paises_dic['Código País'], paises_dic['Nombre País'])) 
    areaf_dic = pd.read_excel(f'{dataDirI}/{path_paises}', sheet_name='Anexo 5', skiprows=8)
    areaf_dic = dict(zip(areaf_dic['Código'], areaf_dic['Descripción'])) 
    comunas_dic = pd.read_csv(f'{dataDirI}/{path_comunas}').rename(columns={'comuna_str':'comuna_paciente_str','Region':'Region_paciente'})


    df = (df
          .fillna({'SEXO': 99})
          .fillna({col: -1 for col in serv_cl_col+areaf_col+columnas_a_convertir}))
    df = (df
          .assign(**{col: df[col].astype(int) for col in serv_cl_col+areaf_col+columnas_a_convertir})
          .pipe(lambda df: df.assign(**{col: pd.Categorical(df[col]) for col in columnas_a_convertir})))
    df = (df
          .assign(sexo_str = lambda df: df['SEXO'].map(sexo),
                  tipo_edad_str = lambda df: df['TIPO_EDAD'].map(tipo_edad_dic),
                  etnia_str = lambda df: df['ETNIA'].map(etnia_dic),
                  pais_origen = lambda df: df['P_ORIGEN'].map(paises_dic),
                  nacional = lambda df: (df['P_ORIGEN']==152).astype(int),
                  COMUNA = lambda df: df['COMUNA'].astype(int)))
    df = (
            df
            .pipe(lambda df: df.assign(**{col: pd.Categorical(df[col]) for col in columnas_a_convertir2}))
            .assign(in_season = lambda df: inseason(df['M_NAC'],df['A_NAC'],df['ANO_ING'])))
    
    df = (df.merge(comunas_dic,how='left',on='COMUNA'))

    df = df.drop(columns=['Unnamed: 0'])

    return df
 
#################################### Función Final de Egresos #####################################
# Función de preparación de la base de datos de egresos
def egresos(path,output_text):
    
    # Cargamos los egresos
    egresos = pre_egresos(path)

    # Filtramos los egresos
    egresos = filtros_egresos(egresos)

    ### Agregamos resultados importantes ###
    # Informción de establecimentos
    egresos = estab(egresos)
    # Información categórica
    egresos = categorical_info(egresos)

    # egresos_full.parquet
    egresos.to_parquet(f'{dataDirO}/egresos_all{output_text}.parquet', engine='pyarrow', compression='snappy')

    #egresos_5.parquet
    df_5 = egresos[egresos['dias_vida']<=365*5+3]
    df_5.to_parquet(f'{dataDirO}/egresos_5{output_text}.parquet', engine='pyarrow', compression='snappy')

    #egresos_respiratory_causes.parquet
    df_respiratory_causes = egresos[egresos['diag_respiratory_causes']==True]
    df_respiratory_causes.to_parquet(f'{dataDirO}/egresos_respiratory_causes{output_text}.parquet', engine='pyarrow', compression='snappy')

# Función final de los datos de egresos con todas las variables necesarias


###################################################################################################
########################################### NACIMIENTOS ###########################################
###################################################################################################

###################################### Funciones auxiliares #######################################
def pre_nacimientos(path):
    nacimientos = pd.read_csv(f'{dataDirI}/{path}',encoding= 'latin1', sep = ';', low_memory=False)
    try:
        nacimientos = (nacimientos
                       .assign(fechaIng=lambda x: pd.to_datetime(x['FECHA_INGRESO'],format='mixed', dayfirst=True, errors='coerce'),
                               fechaNac=lambda x: pd.to_datetime(x['FECHA_NACIMIENTO'],format='mixed', dayfirst=True, errors='coerce'),
                               fechaEgr=lambda x: pd.to_datetime(x['FECHA_EGRESO'],format='mixed', dayfirst=True, errors='coerce'),
                               fechaDef=lambda x: pd.to_datetime(x['FECHA_DEFUNCION'],format='mixed', dayfirst=True, errors='coerce'),
                               fechaInm=lambda x: pd.to_datetime(x['FECHA_INMUNIZACION'],format='mixed', dayfirst=True, errors='coerce'))
                       .assign(DIA_ING=lambda x: x['fechaIng'].dt.day,
                               MES_ING=lambda x: x['fechaIng'].dt.month,
                               ANO_ING=lambda x: x['fechaIng'].dt.year,
                               DIA_EGR=lambda x: x['fechaEgr'].dt.day,
                               MES_EGR=lambda x: x['fechaEgr'].dt.month,
                               ANO_EGR=lambda x: x['fechaEgr'].dt.year,
                               D_NAC=lambda x: x['fechaNac'].dt.day,
                               M_NAC=lambda x: x['fechaNac'].dt.month,
                               A_NAC=lambda x: x['fechaNac'].dt.year))
        return nacimientos
    except:
        
        try:
            nacimientos = (nacimientos
                           .assign(fechaIng=lambda x: pd.to_datetime({'year': x['ANO_ING'], 'month': x['MES_ING'], 'day': x['DIA_ING']},format='%Y-%m-%d'),
                                   fechaNac=lambda x: pd.to_datetime(x['FECHA_NAC'],dayfirst=True),
                                   fechaEgr=lambda x: pd.to_datetime(x['FECHA_EGR'],dayfirst=True),
                                   fechaDef=lambda x: pd.to_datetime(x['FECHA_DEF'],dayfirst=True))
                           .assign(DIA_EGR=lambda x: x['fechaEgr'].dt.day,
                                   MES_EGR=lambda x: x['fechaEgr'].dt.month,
                                   ANO_EGR=lambda x: x['fechaEgr'].dt.year,
                                   D_NAC=lambda x: x['fechaNac'].dt.day,
                                   M_NAC=lambda x: x['fechaNac'].dt.month,
                                   A_NAC=lambda x: x['fechaNac'].dt.year))
            return nacimientos
        except:
            print('nacimientos no tiene el formato conocido')

################################## Función Final de Nacimientos ###################################
# Función final de nacimientos
def nacimientos(path,output_text):
    # Cargamos los nacimeintos
    nacimientos = pre_nacimientos(path)
    # Guardamos la base de datos sin filtrar de nacimientos
    nacimientos.to_parquet(f'{dataDirO}/Nacimientos/nacimientos{output_text}.parquet', engine='pyarrow', compression='snappy')


###################################################################################################
####################################### PACIENTES DE RIESGO #######################################
###################################################################################################

###################################### Funciones auxiliares #######################################

# Función que genera la base de datos de nacimientos con las características de riesgo
def nacimiento_riesgo(path):
    try:
        cols_nac = ['RUN','fechaNac','fechaDef','SEMANAS','PESO',
                'TALLA','D_NAC','M_NAC','A_NAC','RUN_RNI','RUN_M',
                'VACUNADO','MARCA','FECHA_INMUNIZACION','fechaInm']
        nacimientos = (pd
                       .read_parquet(f'{dataDirO}/Nacimientos/{path}')[cols_nac]
                       .sort_values('fechaDef',ascending=True)
                       .drop_duplicates(subset=['RUN']))
    except: 
        try:
            cols_nac = ['RUN','fechaNac','fechaDef','SEMANAS','PESO',
                'TALLA','D_NAC','M_NAC','A_NAC']
            nacimientos = (pd
                           .read_parquet(f'{dataDirO}/Nacimientos/{path}')[cols_nac]
                           .sort_values('fechaDef',ascending=True)
                           .drop_duplicates(subset=['RUN']))
        except:
            nacimientos = pd.read_parquet(f'{dataDirO}/Nacimientos/{path}').drop_duplicates(subset=['RUN'])


    nacimientos = nacimientos.assign(riesgo = lambda x: (x['SEMANAS']<=32)|(x['PESO']<=1500))

    fechas_inicio = {y: pd.Timestamp(y, 4, 1) for y in range(2014, 2025)}
    fechas_fin = {y: pd.Timestamp(y, 10, 1) for y in range(2014, 2025)}

    for y in range(2019, 2025):
        nacimientos[f'catchup2_{y}'] = (
            (nacimientos['fechaNac'] >= fechas_inicio[y - 1]) &
            (nacimientos['fechaNac'] < fechas_fin[y - 1])
        )
        nacimientos[f'catchup_{y}'] = (
            (nacimientos['fechaNac'] >= fechas_fin[y - 1]) &
            (nacimientos['fechaNac'] < fechas_inicio[y])
        )
        nacimientos[f'inseason_{y}'] = (
            (nacimientos['fechaNac'] >= fechas_inicio[y]) &
            (nacimientos['fechaNac'] <= fechas_fin[y])
        )
        nacimientos[f'riesgo_{y}'] = (
            (nacimientos[f'catchup2_{y}'] | nacimientos[f'catchup_{y}'] |nacimientos[f'inseason_{y}']) &
            nacimientos['riesgo']
        )
    return nacimientos

# Función que agrega el resumen de enfermedades a la base de egresos
def enfermedades(df):
    diag_col = [col for col in df.columns if col.startswith("DIAG")]

    # VRS all
    df['vrs0_all'] = df[diag_col].isin(set(rsv_0)).any(axis=1)
    df['vrs1_all'] = df[diag_col].isin(set(rsv_1)).any(axis=1)
    df['vrs2_all'] = df[diag_col].isin(set(rsv_2)).any(axis=1)

    # VRS
    df['vrs0'] = df['DIAG1'].isin(set(rsv_0))
    df['vrs1'] = df['DIAG1'].isin(set(rsv_1))
    df['vrs2'] = df['DIAG1'].isin(set(rsv_2))

    # Respiratory diseases
    diag_codes = pd.read_excel(f'{dataDirI}/IMPACTO/codigos.xlsx', sheet_name='Diag2')
    pneumonia = set(diag_codes['neumonia'].dropna().tolist())
    covid = set(diag_codes['covid'].dropna().tolist())
    influenza = set(diag_codes['influenza'].dropna().tolist())
    acute_bronchitis_bronchiolitis = set(diag_codes['bronquitis'].dropna().tolist())
    obstructive_pulmonary_disease = set(diag_codes['cob'].dropna().tolist())
    upper_respiratory_tract_infection = set(diag_codes['ira_alta'].dropna().tolist())
    other_respiratory_conditions = set(diag_codes['otros'].dropna().tolist())

    df['influenza'] = df[diag_col].isin(influenza).any(axis=1)
    df['pneumonia'] = df[diag_col].isin(pneumonia).any(axis=1)
    df['upper_respiratory_tract_infection'] = df[diag_col].isin(upper_respiratory_tract_infection).any(axis=1)
    df['obstructive_pulmonary_disease'] = df[diag_col].isin(obstructive_pulmonary_disease).any(axis=1)
    df['acute_bronchitis_bronchiolitis'] = df[diag_col].isin(acute_bronchitis_bronchiolitis).any(axis=1)
    df['covid'] = df[diag_col].isin(covid).any(axis=1)
    df['other_respiratory_conditions'] = df[diag_col].isin(other_respiratory_conditions).any(axis=1)
    disease_cols = ['influenza', 'pneumonia',
                'upper_respiratory_tract_infection',
                'obstructive_pulmonary_disease',
                'acute_bronchitis_bronchiolitis',
                'covid',
                'other_respiratory_conditions']
    df['non_respiratory_disease'] = df[disease_cols].sum(axis=1) == 0

    # IRAG
    irag = set(pd.read_csv(f'{dataDirI}/IMPACTO/codigos_IRAG.csv').IRAG.to_list())
    df['irag'] = df[diag_col].isin(irag).any(axis=1)

    # Risk diseases
    df['displasia'] = df[diag_col].isin(set(displacia_pulmonar)).any(axis=1)
    df['cardipatia1'] = df[diag_col].isin(set(cardiopatias)).any(axis=1)
    df['cardipatia2'] = df[diag_col].isin(set(cardiopatias2)).any(axis=1)

    return df

# Funciones que agrega las columnas de riesgo y de elegibilidad a ls bases
def riesgo_elegibilidad_nacimientos(df,fecha_min_nac=(4,1)):
    df = (df
          .assign(bajo_peso = lambda x: (x['PESO'] <= 1500),
                  prematuro = lambda x: (x['SEMANAS'] <= 32))
          .assign(riesgo1=lambda x: (x['SEMANAS'] <= 32) | (x['PESO'] <= 1500) | (x['card1']) | (x['displ']),
                  riesgo2=lambda x: (x['SEMANAS'] <= 32) | (x['PESO'] <= 1500) | (x['card2']) | (x['displ'])))

    fechas_inicio = {y: pd.Timestamp(y, fecha_min_nac[0], fecha_min_nac[1]) for y in range(2018, 2025)}
    fechas_fin = {y: pd.Timestamp(y, 9, 30) for y in range(2019, 2025)}

    for y in range(2019, 2025):

        df[f'riesgo1_{y}'] = (
            (df[f'catchup2_{y}'] | df[f'catchup_{y}'] |df[f'inseason_{y}']) &
            (df['riesgo1'])
        )
        df[f'riesgo2_{y}'] = (
            (df[f'catchup2_{y}'] | df[f'catchup_{y}'] |df[f'inseason_{y}']) &
            (df['riesgo2'])
        )

    return df

def riesgo_elegibilidad_egresos(df,fecha_min_nac=(4,1)):
    fechas_inicio = {y: pd.Timestamp(y, fecha_min_nac[0], fecha_min_nac[1]) for y in range(2018, 2025)}
    fechas_fin = {y: pd.Timestamp(y, 9, 30) for y in range(2019, 2025)}

    df = df.assign(riesgo=lambda x: (x['SEMANAS'] <= 32) | (x['PESO'] <= 1500))
    df = df.assign(riesgo1=lambda x: (x['SEMANAS'] <= 32) | (x['PESO'] <= 1500) | (x['card1']) | (x['displ']),
                   riesgo2=lambda x: (x['SEMANAS'] <= 32) | (x['PESO'] <= 1500) | (x['card2']) | (x['displ']))

    for y in range(2019, 2025):

        df[f'catchup2_{y}'] = (
            (df['fechaNac'] >= fechas_inicio[y - 1]) &
            (df['fechaNac'] < fechas_fin[y - 1])
        )
        df[f'catchup_{y}'] = (
            (df['fechaNac'] >= fechas_fin[y - 1]) &
            (df['fechaNac'] < fechas_inicio[y])
        )
        df[f'inseason_{y}'] = (
            (df['fechaNac'] >= fechas_inicio[y]) &
            (df['fechaNac'] <= fechas_fin[y])
        )
        df[f'riesgo_{y}'] = (
            (df[f'catchup2_{y}'] | df[f'catchup_{y}'] | df[f'inseason_{y}']) &
            df['riesgo']
        )

        df[f'riesgo1_{y}'] = (
            (df[f'catchup2_{y}'] | df[f'catchup_{y}'] |df[f'inseason_{y}']) &
            df['riesgo1']
        )
        df[f'riesgo2_{y}'] = (
            (df[f'catchup2_{y}'] | df[f'catchup_{y}'] |df[f'inseason_{y}']) &
            df['riesgo2']
        )

        df[f'card1_{y}'] = ((df[f'catchup2_{y}'] | df[f'catchup_{y}'] |df[f'inseason_{y}'])&df['card1'])
        df[f'card2_{y}'] = ((df[f'catchup2_{y}'] | df[f'catchup_{y}'] |df[f'inseason_{y}'])&df['card2'])
        df[f'displ_{y}'] = ((df[f'catchup2_{y}'] | df[f'catchup_{y}'] |df[f'inseason_{y}'])&df['displ'])

        df[f'no_eligible_{y}'] = ((df['fechaNac'] >= (df['fechaIng'] - pd.DateOffset(years=4)))
                                  & (df['fechaNac']<fechas_inicio[y-1]))

    return df

##################################### Función Final de Riesgo #####################################

# Función que une la base de egresos con la de nacimientos para detectar los pacientes de riesgo
# Guarda la base de egresos con las caraterísticas de riesgo del paciente y la lista de pacientes de riesgo
def riesgo(path_nacimientos,path_egresos,output_text):

    # Cargamos la base de nacimientos
    nacimientos = nacimiento_riesgo(path_nacimientos)

    # Cargamos la base de egresos
    egresos = pd.read_parquet(f'{dataDirO}/{path_egresos}', engine='pyarrow')

    # Unimos las bases de datos 
    egresos = egresos.merge(nacimientos,on=['RUN','fechaNac'],how='left')

    # Fitros 
    egresos = egresos.drop_duplicates(['RUN','fechaIng'],keep='first')

    # Agregamos las columnas de enfermedades
    egresos = enfermedades(egresos)

    # Obtenemos los pacientes de riesgo
    card1 = egresos.query('cardipatia1').drop_duplicates(subset='RUN')[['RUN','cardipatia1']].rename(columns={'cardipatia1':'card1'})
    card2 = egresos.query('cardipatia2').drop_duplicates(subset='RUN')[['RUN','cardipatia2']].rename(columns={'cardipatia2':'card2'})
    displ = egresos.query('displasia').drop_duplicates(subset='RUN')[['RUN','displasia']].rename(columns={'displasia':'displ'})

    # Agregamos el indicador de riesgo 
    egresos = (egresos
               .merge(card1,on='RUN',how='left')
               .merge(card2,on='RUN',how='left')
               .merge(displ,on='RUN',how='left')
               .fillna({'card1':False,'card2':False,'displ':False}))
    nacimientos = (nacimientos
                   .merge(card1,on='RUN',how='left')
                   .merge(card2,on='RUN',how='left')
                   .merge(displ,on='RUN',how='left')
                   .fillna({'card1':False,'card2':False,'displ':False}))

    nacimientos = riesgo_elegibilidad_nacimientos(nacimientos)
    egresos = riesgo_elegibilidad_egresos(egresos)
    
    # Guardamos las bases de datos 
    egresos.to_parquet(f'{dataDirO}/egresos_riesgo{output_text}.parquet', engine='pyarrow', compression='snappy')
    nacimientos.to_parquet(f'{dataDirO}/Nacimientos/nacimientos_riesgo{output_text}.parquet', engine='pyarrow', compression='snappy')
    # Guardamos la lista de pacientes de riesgo
    nacimientos[['RUN','card1','card2','displ']].drop_duplicates(subset='RUN').to_csv(f'{dataDirO}/lista_pacientes_riesgo{output_text}.csv',index=False)
    



    



    