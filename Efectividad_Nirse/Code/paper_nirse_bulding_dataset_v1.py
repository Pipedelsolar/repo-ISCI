import os
import sys

import pandas as pd
import numpy as np
from datetime import datetime,timedelta

from Data.features import convertir_a_fecha
from Data.features import season
from Data.features import prematuro
from Data.features import retraso
from Data.features import son_similares
from Data.features import mapear_region
from Data.features import obtener_fecha_primer_upc
from Data.data_processing import call_comunas
import config.config as cf
from config import tools


@tools.timeit(return_execution_time=False,execution_time_digits=0)
def call_immunized_discharges(path="/Users/ignasi/Documents/ISCI/Projects/Data/EntregasMinsal/EfectividadNirsevimad/NAC_RNI_EGRESOS_ENTREGA_ISCI_20_08_2024_encr.csv"):
    comunas = call_comunas().rename(columns = {'COMUNA': 'COMUNA_N'}).drop_duplicates(subset='COMUNA_N')
    #check codigo diagnostico
    cod_diag = (
        pd.read_csv('/Users/ignasi/Documents/ISCI/Projects/Data/EntregasMinsal/Codigo_diagnosticos.csv')
        
    )
    IRAG=pd.read_excel('/Users/ignasi/Documents/ISCI/Projects/Data/codigos_IRAG.xlsx',sheet_name='IRAG').IRAG.unique()
    #pd.read_excel('Data/codigos_IRAG.xlsx',sheet_name='Diag2').unique()
    df_weight = pd.read_csv('/Users/ignasi/Documents/ISCI/Projects/nirsevimab/Data/Output/gestation_and_weigh_projected.csv', encoding = "latin-1", sep = ",")

    orden_categoria = ['0-7', '7-14', '14-21', '21-28', '28-42','42-56','56-70','+70', 'no inmunizado']#,'sin categoria'
    orden_categoria_edad = ['menores de 1 mes', '1-2 meses', '2-3 meses', '3-6 meses', '6+ meses']#,'sin información','no ingreso']
    
    return (   
        pd.read_csv(path, sep = ';',encoding = "latin1")#,nrows=1000000
        .rename(columns={
            'FECHA_NACIMIENTO':'FECHA_NAC',
            'FECHA_INGRESO':'FECHA_ING',
            'FECHA_EGRESO': 'FECHA_EGR'
        }
        )
        .merge(comunas,how='left',on ='COMUNA_N')
        .assign(
            # fechaNac = lambda x: pd.to_datetime(x['FECHA_NAC'], format='%d%b%Y'),
            # fechaIng = lambda x: pd.to_datetime(x['FECHA_ING'], format='%d%b%Y'),
            # fechaEgr = lambda x: pd.to_datetime(x['FECHA_EGR'], format='%d%b%Y'),
            # fechaInm = lambda x: pd.to_datetime(x['FECHA_INMUNIZACION'], format='%d%b%Y'),
            fechaNac = lambda x: pd.to_datetime(x['FECHA_NAC'], infer_datetime_format=True),
            fechaIng = lambda x: pd.to_datetime(x['FECHA_ING'], infer_datetime_format=True),
            fechaEgr = lambda x: pd.to_datetime(x['FECHA_EGR'], infer_datetime_format=True),
            fechaInm = lambda x: pd.to_datetime(x['FECHA_INMUNIZACION'], infer_datetime_format=True),
            SEXO =lambda x: x['SEXO'].astype(str).map(cf.sex),
            diag_vrs=lambda df:df[['DIAG1'#, 'DIAG3', 'DIAG4', 'DIAG5', 'DIAG6', 'DIAG7', 'DIAG8', 'DIAG9', 'DIAG10', 'DIAG11'
                        ]].apply(lambda row: row.isin(cf.diagnosticosVRS).any(), axis=1).astype(int).astype(bool),
            diag_respiratory_causes=lambda x: x['DIAG1'].str.startswith('J') | x['DIAG1'].eq('B974') ,
            diag_irag=lambda x: x['DIAG1'].isin(IRAG),
            )
        .pipe(lambda d: d.assign(**{
            col: d[col].apply(lambda x: 1 if x in cf.areasUPC else 0)
            for col in  cf.colsDiagnostico
    
        }))
        .assign(
            estado_inmunizacion=lambda df: np.where( np.logical_and(df['fechaInm'].notna(), ((df['fechaIng'].isna()) | (df['fechaIng'] >= df['fechaInm']))),
                                                     'inmunizado',
                                                     'not inmunizado'),
            NOMBRE_REGION = lambda df: df['NOMBRE_REGION'].str.replace('Región', '').str.strip()
                                            .fillna('Ignorada')
                                            .apply(lambda x: mapear_region(x, cf.dict_regiones)),
            fechaIng_week = lambda df: df.fechaIng.dt.isocalendar().week,
        )
        .assign(
            #update Inm
            fechaInm_cox=lambda df: df.apply(lambda row: pd.NaT if ((row['fechaIng'] - row['fechaInm']).days <= 7 and row['diag_vrs'] == 1 and not pd.isna(row['fechaInm'])) else row['fechaInm'], axis=1),
            #check vrs
            fechaIng_vrs= lambda df: df.apply(lambda row: pd.NaT if (not row['diag_vrs'] and not pd.isna(row['fechaIng'])) else row['fechaIng'], axis=1),

            Macrozona=lambda df: df['NOMBRE_REGION'].map(cf.region_a_macrozona),
            Macrozona2=lambda df: df['NOMBRE_REGION'].map(cf.region_a_macrozona2),
            group = lambda df: np.where(df['fechaNac']< pd.to_datetime("2024-04-01"), "CATCH_UP", "SEASONAL"),
            #season = lambda x: x['fechaNac'].apply(season),
            #prematuro = lambda x: x['SEMANAS'].apply(prematuro),
        )
        .assign(
             dias_Inm=lambda x:  np.where(
                                        x['fechaIng'].isna() | x['fechaInm'].isna(), 
                                         np.nan,
                                        (x['fechaIng'] - x['fechaInm']).dt.days
                                    ),
             edad_meses_Ing=lambda x: np.where(
                                        x['fechaIng'].isna() | x['fechaNac'].isna(), 
                                         np.nan,
                                        ((x['fechaIng'] - x['fechaNac']).dt.days // 30)
                                    )
        )
        .assign(
        categoria_dias_Inm=lambda x:  pd.Categorical(np.select(
            [
                #x['fechaInm'].isna(),  # No tiene fecha de inmunización
                (x['dias_Inm'] >= 0) & (x['dias_Inm'] <= 7),
                (x['dias_Inm'] > 7) & (x['dias_Inm'] <= 14),
                (x['dias_Inm'] > 14) & (x['dias_Inm'] <= 21),
                (x['dias_Inm'] > 21) & (x['dias_Inm'] <= 28),
                (x['dias_Inm'] > 28) & (x['dias_Inm'] <= 42),
                (x['dias_Inm'] > 42) & (x['dias_Inm'] <= 56),
                (x['dias_Inm'] > 56) & (x['dias_Inm'] <= 70),

                (x['dias_Inm'] > 70)
            ],
            [
                '0-7',
                '7-14',
                '14-21',
                '21-28',
                '28-42',
                '42-56',
                '56-70',
                '+70'
            ],
            default= 'no inmunizado'
            ), categories=orden_categoria, ordered=True),
        categoria_edad_meses_Ing=lambda x: pd.Categorical( np.select(
            [   
                x['edad_meses_Ing'].isna(),
                x['edad_meses_Ing'] < 1,
                (x['edad_meses_Ing'] >= 1) & (x['edad_meses_Ing'] < 2),
                (x['edad_meses_Ing'] >= 2) & (x['edad_meses_Ing'] < 3),
                (x['edad_meses_Ing'] >= 3) & (x['edad_meses_Ing'] < 6),
                x['edad_meses_Ing'] >= 6
            ],
            [   
                'no ingreso',
                'menores de 1 mes',
                '1-2 meses',
                '2-3 meses',
                '3-6 meses',
                '6+ meses'
            ],
            default='sin información'
            ), categories=orden_categoria_edad, ordered=True)
        )
    )