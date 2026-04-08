import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import sys
import os
from dateutil.relativedelta import relativedelta
parent_dir = os.getcwd()
if parent_dir not in sys.path: sys.path.insert(0,parent_dir+'/src')
from config import config as cf
from config import functions as f
folder_path_input = cf.dataDirI2
folder_path = cf.dataDirO2
tribu_priv = (pd.read_excel(cf.dataDirI2+cf.tribu_priv,skiprows=2))
tribu_priv = (
    tribu_priv
    .assign(origen_estab=lambda df: 'private')
    .rename(columns={'SEREMI de Salud':'NombreSeremiSalud','Codigo Establecimiento':'ESTAB','Nombre Establecimiento': 'NombreEstablecimiento'})
    .assign(NombreSeremiSalud = lambda df: df['NombreSeremiSalud'].ffill())
    .dropna(subset= ['ESTAB'])
    .assign(ESTAB=lambda df: df['ESTAB'].astype('int'))
    .assign(NombreSeremiSalud = lambda df: df['NombreSeremiSalud'].str.replace('SEREMI ','',regex=False))
    .assign(NombreSeremiSalud = lambda df: df['NombreSeremiSalud'].str.replace('de ','',regex=False))
    .assign(NombreSeremiSalud = lambda df: df['NombreSeremiSalud'].str.strip())
    .assign(Region_SS = lambda df: df['NombreSeremiSalud'].map(cf.trib_servicios_dic).map(cf.trib_servicios_reg_dic)))

tribu_pub = (pd.read_excel(cf.dataDirI2+cf.tribu_pub,skiprows=2))
tribu_pub = (
    tribu_pub
    .assign(origen_estab=lambda df: 'public')
    .rename(columns={'Servicio de Salud': 'NombreServicioSalud', 'Codigo Establecimiento':'ESTAB','Nombre Establecimiento': 'NombreEstablecimiento'})
    .assign(NombreServicioSalud = lambda df: df['NombreServicioSalud'].ffill())
    .assign(NombreServicioSalud = lambda df: df['NombreServicioSalud'].map(cf.trib_servicios_dic))
    .assign(Region_SS = lambda trib: trib['NombreServicioSalud'].map(cf.trib_servicios_reg_dic))
    .dropna(subset= ['ESTAB'])
    .assign(ESTAB=lambda df: df['ESTAB'].astype('int'))
    )

trib = (pd.concat([tribu_priv,tribu_pub]))
trib = (
    trib
    .assign(trib_suf = lambda trib:  trib['Jun'] <= 0.8 * trib['May'])
    .assign(ESTAB = lambda trib:  trib['ESTAB'].astype('category'))
    .assign(Macrozona1_SS = lambda trib: trib['Region'].map(cf.region_a_macrozona_V2))
    .assign(Macrozona2_SS = lambda trib: trib['Region'].map(cf.region_a_macrozona2_V2))
    .assign(trib_rech = lambda trib: trib['ESTAB'].isin(cf.filter_estab))
)
trib.to_csv(f'{folder_path}/trib.csv')
establecimientos = (pd
                     .read_excel(folder_path_input+cf.establecimientos)
                     .rename(columns={'Código Comuna':'COMUNA','Nombre Comuna':'comuna_str', 'Nombre Región': 'Region', 'COD_VIG':'ESTAB'})
                     .dropna(subset={'ESTAB'})
                     .assign(ESTAB = lambda establecimientos: establecimientos['ESTAB'].astype(int)))

comunas_dic = (establecimientos[['COMUNA','comuna_str','Region']]
               .assign(Region = lambda establecimientos: establecimientos['Region'].map(cf.dict_regiones2))
               .drop_duplicates())

comunas_dic.to_csv(f'{folder_path}/comunas_dic.csv')

estab_loc = (establecimientos[['ESTAB','COMUNA','comuna_str','Region']]
             .rename(columns={'COMUNA':'COMUNA_ESTAB','comuna_str':'comuna_estab_str','Region':'Region_ESTAB'})
             .drop_duplicates()
             .assign(Region = lambda establecimientos: establecimientos['Region_ESTAB'].map(cf.dict_regiones2))
             )
estab_loc.to_csv(f'{folder_path}/estab_loc.csv')
egresos_all = pd.read_csv(cf.dataDirI2+cf.egresos,encoding= 'latin1', sep = '|', low_memory=False)

trib = pd.read_csv(f'{folder_path}/trib.csv')
paises_dic = pd.read_excel(cf.dataDirI2+cf.hospitalaria, sheet_name='Anexo 3', skiprows=6)
paises_dic = dict(zip(paises_dic['Código País'], paises_dic['Nombre País'])) 
areaf_dic = pd.read_excel(cf.dataDirI2+cf.hospitalaria, sheet_name='Anexo 5', skiprows=8)
areaf_dic = dict(zip(areaf_dic['Código'], areaf_dic['Descripción'])) 
serv_cl_col = list([f'SERC_{x}_TRAS' for x in range(1,10)])+['SER_CLIN_I','SERC_EGR']
areaf_col = list([f'AREAF_{x}_TRAS' for x in range(1,10)])+['AREAF_1_TRAS','AREAF_EGR']
columnas_a_convertir = ['ESTAB',
                         'ServicioSalud',
                         'Seremi',
                         'SEXO',
                         'TIPO_EDAD',
                         'ETNIA',
                         'P_ORIGEN',
                         'COMUNA',
                         'PREVI',
                         'COND_EGR']
columnas_a_convertir2 = ['sexo_str','tipo_edad_str','etnia_str','pais_origen','nacional','COMUNA']
estab_loc = (pd
             .read_csv(f'{folder_path}/estab_loc.csv')[['ESTAB','COMUNA_ESTAB','comuna_estab_str','Region']]
             .rename(columns={'Region':'Region_ESTAB'}))
df = (
        egresos_all
        .rename(columns={'RUT':'RUN'})
        .assign(
                fechaNac =  lambda df: pd.to_datetime({'year': df['A_NAC'], 'month': df['M_NAC'], 'day': df['D_NAC']}, format='%Y-%m-%d'),
                fechaIng = lambda df: pd.to_datetime({'year': df['ANO_ING'], 'month': df['MES_ING'], 'day': df['DIA_ING']}, format='%Y-%m-%d')
        )
)
df = df.assign(dias_vida =lambda x: ((x['fechaIng'] - x['fechaNac']).dt.days))
df = (
        df
        .assign( 
                fechaEgr = lambda df: pd.to_datetime({'year': df['ANO_EGR'], 'month': df['MES_EGR'], 'day': df['DIA_EGR']}, format='%Y-%m-%d'),
                duracion = lambda df: (df['fechaEgr']-df['fechaIng']).dt.days.apply(lambda x: max(1, x)),
                epiweek=lambda x: x['fechaIng'].dt.isocalendar().week,
                year=lambda x: x['fechaIng'].dt.isocalendar().year,
                first_letter =  lambda x: x.DIAG1.str[0],
                diag_vrs=lambda df:df[['DIAG1'
                            ]].apply(lambda row: row.isin(cf.diagnosticosVRS).any(), axis=1).astype(int).astype(bool),
                diag_vrs_diag3=lambda df:df[['DIAG3'
                            ]].apply(lambda row: row.isin(cf.diagnosticosVRS).any(), axis=1).astype(int).astype(bool),
                diag_respiratory_causes=lambda x: x['DIAG1'].str.startswith('J') | x['DIAG1'].eq('B974')
        ))
df.to_parquet(folder_path+'/egresos_sin_filtro.parquet', engine='pyarrow', compression='snappy')
df = pd.read_parquet(folder_path+'/egresos_sin_filtro.parquet')
df = df[(df['year'] >= 2018) & 
        (~df['ESTAB'].isna()) & 
        (df['dias_vida'] >= 0) & 
        (df['DIAS_ESTAD'] >= 0) &
        (~df['fechaIng'].isna())]
df = (df
      .merge(trib[['ESTAB',
                   'NombreEstablecimiento',
                   'origen_estab',
                   'NombreServicioSalud',
                   'NombreSeremiSalud',
                   'Region_SS',
                   'Macrozona1_SS',
                   'Macrozona2_SS',
                   'trib_rech']],
             how='left',
             on='ESTAB')
      .fillna({'trib_rech': True}))

df = (df.merge(estab_loc,how='left',on='ESTAB'))
df = (
        df
        .fillna({'SEXO': 99})
        .fillna({col: -1 for col in serv_cl_col+areaf_col+columnas_a_convertir}))
df = (
        df
        .assign(**{col: df[col].astype(int) for col in serv_cl_col+areaf_col+columnas_a_convertir})
        .pipe(lambda df: df.assign(**{col: pd.Categorical(df[col]) for col in columnas_a_convertir})))
df = (
        df
        .assign(sexo_str = lambda df: df['SEXO'].map(cf.sexo),
                tipo_edad_str = lambda df: df['TIPO_EDAD'].map(cf.tipo_edad_dic),
                etnia_str = lambda df: df['ETNIA'].map(cf.etnia_dic),
                pais_origen = lambda df: df['P_ORIGEN'].map(paises_dic),
                nacional = lambda df: (df['P_ORIGEN']==152).astype(int),
                COMUNA = lambda df: df['COMUNA'].astype(int)))
df = (
        df
        .pipe(lambda df: df.assign(**{col: pd.Categorical(df[col]) for col in columnas_a_convertir2}))
        .assign(in_season = lambda df: f.inseason(df['M_NAC'],df['A_NAC'],df['ANO_ING'])))

comunas_dic = pd.read_csv(f'{folder_path}/comunas_dic.csv').rename(columns={'comuna_str':'comuna_paciente_str','Region':'Region_paciente'})
df = (df.merge(comunas_dic,how='left',on='COMUNA'))
df = df.drop(columns=['Unnamed: 0'])
# egresos_full.parquet
df.to_parquet(folder_path+'/egresos_all.parquet', engine='pyarrow', compression='snappy')

#egresos_5.parquet
df_5 = df[df['dias_vida']<=365*5+3]
df_5.to_parquet(folder_path+'/egresos_5.parquet', engine='pyarrow', compression='snappy')
#egresos_respiratory_causes.parquet
df_respiratory_causes = df[df['diag_respiratory_causes']==True]
df_respiratory_causes.to_parquet(f'{folder_path}/egresos_respiratory_causes.parquet', engine='pyarrow', compression='snappy')
egresos_all = pd.read_csv(cf.dataDirI2+cf.egresos2024,encoding= 'latin1', sep = '|', low_memory=False)

df = (
        egresos_all
        .rename(columns={'RUT':'RUN'})
        .assign(
                fechaNac =  lambda df: pd.to_datetime({'year': df['A_NAC'], 'month': df['M_NAC'], 'day': df['D_NAC']}, format='%Y-%m-%d'),
                fechaIng = lambda df: pd.to_datetime({'year': df['ANO_ING'], 'month': df['MES_ING'], 'day': df['DIA_ING']}, format='%Y-%m-%d')
        )
)
df = df.assign(dias_vida =lambda x: ((x['fechaIng'] - x['fechaNac']).dt.days))
df = (
        df
        .assign( 
                fechaEgr = lambda df: pd.to_datetime({'year': df['ANO_EGR'], 'month': df['MES_EGR'], 'day': df['DIA_EGR']}, format='%Y-%m-%d'),
                duracion = lambda df: (df['fechaEgr']-df['fechaIng']).dt.days.apply(lambda x: max(1, x)),
                epiweek=lambda x: x['fechaIng'].dt.isocalendar().week,
                year=lambda x: x['fechaIng'].dt.isocalendar().year,
                first_letter =  lambda x: x.DIAG1.str[0],
                diag_vrs=lambda df:df[['DIAG1'
                            ]].apply(lambda row: row.isin(cf.diagnosticosVRS).any(), axis=1).astype(int).astype(bool),
                diag_vrs_diag3=lambda df:df[['DIAG3'
                            ]].apply(lambda row: row.isin(cf.diagnosticosVRS).any(), axis=1).astype(int).astype(bool),
                diag_respiratory_causes=lambda x: x['DIAG1'].str.startswith('J') | x['DIAG1'].eq('B974')
        ))
df.to_parquet(folder_path+'/egresos_sin_filtro_2024.parquet', engine='pyarrow', compression='snappy')
df = pd.read_parquet(folder_path+'/egresos_sin_filtro_2024.parquet')
df = df[(df['year'] >= 2018) & 
        (~df['ESTAB'].isna()) & 
        (df['dias_vida'] >= 0) & 
        (df['DIAS_ESTAD'] >= 0) &
        (~df['fechaIng'].isna())]
trib = pd.read_csv(f'{folder_path}/trib.csv')
paises_dic = pd.read_excel(cf.dataDirI2+cf.hospitalaria, sheet_name='Anexo 3', skiprows=6)
paises_dic = dict(zip(paises_dic['Código País'], paises_dic['Nombre País'])) 
areaf_dic = pd.read_excel(cf.dataDirI2+cf.hospitalaria, sheet_name='Anexo 5', skiprows=8)
areaf_dic = dict(zip(areaf_dic['Código'], areaf_dic['Descripción'])) 
serv_cl_col = list([f'SERC_{x}_TRAS' for x in range(1,10)])+['SER_CLIN_I','SERC_EGR']
areaf_col = list([f'AREAF_{x}_TRAS' for x in range(1,10)])+['AREAF_1_TRAS','AREAF_EGR']
columnas_a_convertir = ['ESTAB',
                         'ServicioSalud',
                         'Seremi',
                         'SEXO',
                         'TIPO_EDAD',
                         'ETNIA',
                         'P_ORIGEN',
                         'COMUNA',
                         'PREVI',
                         'COND_EGR']
columnas_a_convertir2 = ['sexo_str','tipo_edad_str','etnia_str','pais_origen','nacional','COMUNA']
estab_loc = (pd
             .read_csv(f'{folder_path}/estab_loc.csv')[['ESTAB','COMUNA_ESTAB','comuna_estab_str','Region']]
             .rename(columns={'Region':'Region_ESTAB'}))
df = (df
      .merge(trib[['ESTAB',
                   'NombreEstablecimiento',
                   'origen_estab',
                   'NombreServicioSalud',
                   'NombreSeremiSalud',
                   'Region_SS',
                   'Macrozona1_SS',
                   'Macrozona2_SS',
                   'trib_rech']],
             how='left',
             on='ESTAB')
      .fillna({'trib_rech': True}))

df = (df.merge(estab_loc,how='left',on='ESTAB'))
df = (
        df
        .fillna({'SEXO': 99})
        .fillna({col: -1 for col in serv_cl_col+areaf_col+columnas_a_convertir}))
df = (
        df
        .assign(**{col: df[col].astype(int) for col in serv_cl_col+areaf_col+columnas_a_convertir})
        .pipe(lambda df: df.assign(**{col: pd.Categorical(df[col]) for col in columnas_a_convertir})))
df = (
        df
        .assign(sexo_str = lambda df: df['SEXO'].map(cf.sexo),
                tipo_edad_str = lambda df: df['TIPO_EDAD'].map(cf.tipo_edad_dic),
                etnia_str = lambda df: df['ETNIA'].map(cf.etnia_dic),
                pais_origen = lambda df: df['P_ORIGEN'].map(paises_dic),
                nacional = lambda df: (df['P_ORIGEN']==152).astype(int),
                COMUNA = lambda df: df['COMUNA'].astype(int)))
df = (
        df
        .pipe(lambda df: df.assign(**{col: pd.Categorical(df[col]) for col in columnas_a_convertir2}))
        .assign(in_season = lambda df: f.inseason(df['M_NAC'],df['A_NAC'],df['ANO_ING'])))

comunas_dic = pd.read_csv(f'{folder_path}/comunas_dic.csv').rename(columns={'comuna_str':'comuna_paciente_str','Region':'Region_paciente'})
df = (df.merge(comunas_dic,how='left',on='COMUNA'))
df = df.drop(columns=['Unnamed: 0'])
# egresos_full.parquet
df.to_parquet(folder_path+'/egresos_all_2024.parquet', engine='pyarrow', compression='snappy')

#egresos_5.parquet
df_5 = df[df['dias_vida']<=365*5+3]
df_5.to_parquet(folder_path+'/egresos_5_2024.parquet', engine='pyarrow', compression='snappy')
#egresos_respiratory_causes.parquet
df_respiratory_causes = df[df['diag_respiratory_causes']==True]
df_respiratory_causes.to_parquet(f'{folder_path}/egresos_respiratory_causes_2024.parquet', engine='pyarrow', compression='snappy')

