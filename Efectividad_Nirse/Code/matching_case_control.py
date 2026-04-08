import os 
import sys
from pathlib import Path
import pickle
import sklearn
from difflib import SequenceMatcher
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from scipy.spatial.distance import mahalanobis
from IPython.display import display
import pandas as pd
import numpy as np
from scipy import stats
# from preproces_prod4_2025_update import pre_filtred, dias, meses, filtros_IH_new
# from preproces_prod3 import *
from pulp import LpMaximize, LpProblem, LpVariable, lpSum
from gurobipy import Model, GRB
import time 
    
local_path = os.getcwd()
code_root = os.path.abspath(os.path.join(local_path, '../..', 'inv'))

if code_root not in sys.path:
    sys.path.insert(0, code_root)

from statsmodels.discrete.conditional_models import ConditionalLogit

np.random.seed(42)
path_actual = Path.cwd()
path_data = path_actual.parent / 'Data'

dias = list(range(7, (40-13)*7, 7)) 
meses = [1,2,3,4,5,6]
diagnosticosVRS = ['J121', 'J205','J219', 'J210', 'B974']##['J120','J122','J123','J204','J206','J207','J211'] #['J121', 'J205','J219', 'J210', 'B974'] #['A099', 'A080', 'A083', 'A082', 'A084', 'A090', 'A081', 'A085'] pachecos_code, leos_code ['J120','J122','J123','J204','J206','J207','J211'] 
diagnosticos_upc = [406, 412, 415, 405, 411, 414, 310, 311, 312, 320, 323, 324]
months_list = ['Abril', 'Mayo', 'Junio', 'Julio', 'Agosto', 'Septiembre', 'Octubre']

reg_print=[]

cache = {'DESCONOCIDO': None,
 'Región De Antofagasta': 'ANTOFAGASTA',
 'Región De Arica Parinacota': 'ARICA Y PARINACOTA',
 'Región De Atacama': 'ATACAMA',
 'Región De Aysén del General Carlos Ibañez del Campo': 'AISEN',
 'Región De Coquimbo': 'COQUIMBO',
 'Región De La Araucanía': 'ARAUCANIA',
 'Región De Los Lagos': 'LOS LAGOS',
 'Región De Los Ríos': 'LOS RIOS',
 'Región De Magallanes y de la Antártica Chilena': 'MAGALLANES Y ANTARTICA',
 'Región De Tarapacá': 'TARAPACA',
 'Región De Valparaíso': 'VALPARAISO',
 'Región De Ñuble': 'NUBLE',
 'Región Del Bíobío': 'BIOBIO',
 "Región Del Libertador Gral. B. O'Higgins": "O'HIGGINS",
 'Región Del Maule': 'MAULE',
 'Región Metropolitana de Santiago': 'METROPOLITANA'}

def mapear_region(region, regiones_dict):
    if region in cache:  # Revisar si ya está en caché
        reg_print.append(f"{region}---->{cache[region]}")
        return cache[region]

    mejor_similitud = 0
    mejor_region = None
    for key, value in regiones_dict.items():
        similitud = SequenceMatcher(None, region, key).ratio()
        if similitud > mejor_similitud:
            mejor_similitud = similitud
            mejor_region = value

    resultado = mejor_region if mejor_similitud > 0 else None
    cache[region] = resultado
    reg_print.append(f"{region}---->{resultado}")
    return resultado

region_a_macrozona2 = {
    "ARICA Y PARINACOTA": "Norte",
    "TARAPACA": "Norte",
    "ANTOFAGASTA": "Norte",
    "ATACAMA": "Norte",  
    "COQUIMBO": "Centro",  
    "VALPARAISO": "Centro",
    "METROPOLITANA": "Centro",
    "O'HIGGINS": "Centro",
    "MAULE": "Centro",  
    "NUBLE": "Centro",
    "BIOBIO": "Centro",
    "ARAUCANIA": "Sur", 
    "LOS RIOS": "Sur",
    "LOS LAGOS": "Sur",
    "AISEN": "Sur",
    "MAGALLANES Y ANTARTICA": "Sur"
} 

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
'DESCONOCIDO' : None,
'De Atacama': 'ATACAMA'
}

def pre_filtred(df_name, lrti_name ='LRTI_Flag'):
    
    path_actual = Path.cwd()
    path_data = path_actual.parent / 'Data'
    
    #Lectura bases
    comunas = (pd.read_excel(path_data/"comunas.xlsx")
        .rename(columns = {'C_COM': 'COMUNA_N','NOM_REG':'NOMBRE_REGION'})
        #.assign(COMUNA = lambda x: x.COMUNA_N)
        .drop_duplicates(subset='COMUNA_N')
    )
    
    df_urba_rural_percent = pd.read_csv(path_data/'df_urba_rural_percent.csv')
    df_com_poor = pd.read_csv(path_data/'df_com_poor.csv')
    dfpeso = pd.read_csv(path_data / 'gestation_and_weigh_projected.csv', encoding = "latin-1", sep = ",")
    
    df_initial = pd.read_csv(path_data / df_name, encoding = "latin-1", sep = ";")
    print(f'n_rows_inicial= {df_initial.shape[0]}')
    
    #Merging
    df = (df_initial
          .copy()
          .merge(comunas,how='left',on ='COMUNA_N') # COMUNA_N
          .merge(df_urba_rural_percent, on='COMUNA', how='left')
          .merge(df_com_poor,on='COMUNA',how='left')
          .merge(dfpeso[['SEMANAS', 'p_00001_lognormal', 'p_99999_lognormal']], on='SEMANAS',how='left')
    )
    
    cols_diagnostico = ['AREA_FUNC_I','AREAF_1_TRAS', 'AREAF_2_TRAS', 'AREAF_3_TRAS', 'AREAF_4_TRAS', 'AREAF_5_TRAS', 'AREAF_6_TRAS', 'AREAF_7_TRAS', 'AREAF_8_TRAS', 'AREAF_9_TRAS','AREAF_EGR']

    tras_date = {'AREA_FUNC_I': 'fechaIng_any','AREAF_1_TRAS':'fecha_tras_1', 'AREAF_2_TRAS':'fecha_tras_2', 'AREAF_3_TRAS':'fecha_tras_3', 'AREAF_4_TRAS':'fecha_tras_4', 'AREAF_5_TRAS':'fecha_tras_5'
                , 'AREAF_6_TRAS':'fecha_tras_6', 'AREAF_7_TRAS':'fecha_tras_7', 'AREAF_8_TRAS':'fecha_tras_8', 'AREAF_9_TRAS':'fecha_tras_9','AREAF_EGR': 'fechaEgr'}

    for i in range(1, 10):
        year_col = f'ANO_{i}_TRAS'
        month_col = f'MES_{i}_TRAS'
        day_col = f'DIA_{i}_TRAS'
        date_col = f'fecha_tras_{i}'
        df[date_col] = pd.to_datetime({'year': df[year_col], 'month': df[month_col], 'day': df[day_col]}, format='%Y-%m-%d')
        
    for col in cols_diagnostico:
        df[col] = df[col].apply(lambda x: 1 if x in diagnosticos_upc else 0)

    def obtener_fecha_primer_upc(row):
        for col in cols_diagnostico:
            if row[col] == 1:
                fecha_col = tras_date[col]
                return row[fecha_col]
        return pd.NaT
    
    IRAG = ['J09X', 'J100', 'J101', 'J108', 'J110', 'J111', 'J118', 'J120',
       'J121', 'J122', 'J123', 'J128', 'J129', 'J13X', 'J14X', 'J150',
       'J151', 'J152', 'J153', 'J154', 'J155', 'J156', 'J157', 'J158',
       'J159', 'J160', 'J168', 'J180', 'J181', 'J182', 'J188', 'J189',
       'J200', 'J201', 'J202', 'J203', 'J204', 'J205', 'J206', 'J207',
       'J208', 'J209', 'J210', 'J211', 'J218', 'J219', 'J22X', 'J40X',
       'J410', 'J411', 'J418', 'J42X', 'J430', 'J431', 'J432', 'J438',
       'J439', 'J440', 'J441', 'J448', 'J449', 'J450', 'J451', 'J458',
       'J459', 'J46X', 'J47X', 'J60X', 'J61X', 'J620', 'J628', 'J630',
       'J631', 'J632', 'J633', 'J634', 'J635', 'J638', 'J64X', 'J65X',
       'J660', 'J661', 'J662', 'J668', 'J670', 'J671', 'J672', 'J673',
       'J674', 'J675', 'J676', 'J677', 'J678', 'J679', 'J680', 'J681',
       'J682', 'J683', 'J684', 'J688', 'J689', 'J690', 'J691', 'J698',
       'J700', 'J701', 'J702', 'J703', 'J704', 'J708', 'J709', 'J80X',
       'J81X', 'J82X', 'J840', 'J841', 'J848', 'J849', 'J850', 'J851',
       'J852', 'J853', 'J860', 'J869', 'J90X', 'J920', 'J929', 'J930',
       'J931', 'J938', 'J939', 'J940', 'J941', 'J942', 'J948', 'J949',
       'J960', 'J961', 'J969', 'J980', 'J981', 'J982', 'J983', 'J984',
       'J985', 'J986', 'J988', 'J989']  #'U071', 'U072'    
    
    ira_alta = ["J00X", "J010", "J011", "J012", "J013", "J014", "J018", "J019", 
                "J020", "J028", "J029", "J030", "J038", "J039", "J040", "J041", 
                "J042", "J050", "J051", "J060", "J068", "J069"]
    
    def filter_lrti_codes(df, cols):
        def is_in_range(code):
            if pd.isnull(code):  # Si el valor es NaN o nulo, retorna False
                return False
            if code in ('J22X','J121', 'J205', 'J210', 'J219', 'B974'):
                return True
            if code.startswith('J') and len(code) >= 4 and code[1:4].isdigit():  # Verificamos si hay números después de 'J'
                return 209 <= int(code[1:4]) <= 229 

            else: 
                return False
         
        return df[cols].applymap(is_in_range).any(axis=1)
    
    df = (df
          .assign(
              #FECHAS: format='%d%b%Y', infer_datetime_format=True)
              
              fecha_nac = lambda x: pd.to_datetime(x['FECHA_NACIMIENTO'], format='mixed'),
              fechaIng_any = lambda x: pd.to_datetime(x['FECHA_INGRESO'], format='mixed'),
              fechaEgr = lambda x: pd.to_datetime(x['FECHA_EGRESO'], format='mixed'),
              FECHA_INMUNIZACION = lambda x: pd.to_datetime(x['FECHA_INMUNIZACION'], format='mixed'),
              fechaInm = lambda x: pd.to_datetime(x['FECHA_INMUNIZACION'], format='mixed'), #, infer_datetime_format=True),
              
              #VRS
              VRS_D1 = lambda x: np.where(x[['DIAG1']].isin(diagnosticosVRS).any(axis=1), 1, 0),
              VRS_D1y3 = lambda x: np.where(x[['DIAG1','DIAG3']].isin(diagnosticosVRS).any(axis=1), 1, 0),
              VRS_Dall = lambda x: np.where(x[['DIAG1','DIAG2','DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9',
                                               'DIAG10','DIAG11']].isin(diagnosticosVRS).any(axis=1), 1, 0),
              
              diag_irag = lambda x: x['DIAG1'].isin(IRAG), 
              diag_ira_alta = lambda x: x['DIAG1'].isin(ira_alta), 
              LRTI_Flag = lambda x: filter_lrti_codes(x, ['DIAG1']),
              LRTI_all_j = lambda x: x['DIAG1'].str.startswith('J'),
              
              LRTI_Flag_Dall = lambda x: filter_lrti_codes(x, ['DIAG1','DIAG2','DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9',
                                               'DIAG10','DIAG11']),
              
              fechaIng_vrs = lambda x: x.apply(lambda row: row['fechaIng_any'] if row['VRS_D1'] == 1 else pd.NaT, axis=1),
              fechaIng_LRTI = lambda x: x.apply(lambda row: row['fechaIng_any'] if row[lrti_name] == True else pd.NaT, axis=1),          
              fechaIng_vrs_Dall = lambda x: x.apply(lambda row: row['fechaIng_any'] if row['VRS_Dall'] == 1 else pd.NaT, axis=1),
              fechaIng_LRTI_Dall = lambda x: x.apply(lambda row: row['fechaIng_any'] if row['LRTI_Flag_Dall'] == True else pd.NaT, axis=1),
              
              #################################################################################################################################################
            #   fechaIng_vrs = lambda x: x.fechaIng_vrs_Dall,
            #   fechaIng_LRTI = lambda x: x.fechaIng_LRTI_Dall,
              #################################################################################################################################################
            
              eleg_group = lambda x: np.where(x.fecha_nac.dt.month.between(4, 9), 'seasonal' + '_' + (x.fecha_nac.dt.year).astype(str),
                                                np.where((x.fecha_nac.dt.month == 3) & (x.fecha_nac.dt.year == 2025), 'seasonal_2025',
                                                         np.where((x.fecha_nac.dt.month >= 10) & (x.fecha_nac.dt.year == 2024), 'exceso_seasonal_2024_and_catchup_2025',
                                                                  np.where(x.fecha_nac.dt.month <= 3, 'catchup' + '_' + (x.fecha_nac.dt.year).astype(str),
                                                                           'catchup' + '_' + (x.fecha_nac.dt.year + 1).astype(str)
                                                                           )
                                                                  )
                                                         )
                                                ),
              eleg_2025 = lambda x: np.where((x.fecha_nac.dt.month.between(3, 9)) & (x.fecha_nac.dt.year == 2025), 'SEASONAL',
                                             np.where((x.fecha_nac.dt.month <= 2) & (x.fecha_nac.dt.year == 2025), 'CATCH_UP',
                                                      np.where((x.fecha_nac.dt.month >= 10 ) & (x.fecha_nac.dt.year == 2024), 'CATCH_UP',
                                                               'no_elegible_2025'
                                                               )
                                                      )
                                             ),
              eleg_2024 = lambda x: np.where((x.fecha_nac.dt.month.between(4, 9)) & (x.fecha_nac.dt.year == 2024), 'SEASONAL',
                                             np.where((x.fecha_nac.dt.month <= 3) & (x.fecha_nac.dt.year == 2024), 'CATCH_UP',
                                                      np.where((x.fecha_nac.dt.month >= 10 ) & (x.fecha_nac.dt.year == 2023), 'CATCH_UP',
                                                               'no_elegible_2024'
                                                               )
                                                      )
                                             ),
              group = lambda x: np.where(x.fecha_nac.dt.month.between(4, 9), "SEASONAL", "CATCH_UP"),
              sexo = lambda x: np.where(x['SEXO']==1.0,1,0),
              prematuro_extremo = lambda x: np.where((x['SEMANAS']>=22) & (x['SEMANAS']<=27),1,0),
              muy_prematuro = lambda x: np.where((x['SEMANAS']>=28) & (x['SEMANAS']<=32),1,0),
              prematuro_moderado = lambda x: np.where((x['SEMANAS']>=33) & (x['SEMANAS']<=36),1,0),
              prematuro = lambda x: np.where((x['SEMANAS']<=35),1,0),  ### CAMBIE ESTA WEA 15-05-2025 de 36 a 35
              atypic_mom_age = lambda x: np.where((x.EDAD_M>=45) | (x.EDAD_M<=19), 1, 0),
              region = lambda x: x['NOMBRE_REGION'].fillna('DESCONOCIDO').apply(lambda x: mapear_region(x, regiones) if isinstance(x, str) else x),
              Macrozona2 = lambda x: x['region'].map(region_a_macrozona2),
              
              mes_nac_name = lambda x: x.fecha_nac.dt.month_name(),  
              semana_nac = lambda x: x.fecha_nac.dt.isocalendar().week.astype(int),
              year_nac = lambda x: x.fecha_nac.dt.year.astype(int),

              
              #UPC things
              cama = lambda x: np.where(x[cols_diagnostico].eq(1).any(axis=1),'UPC', ""),
              fecha_upc = lambda x: x.apply(obtener_fecha_primer_upc, axis=1),##################FECHA UPC##############
              fecha_upc_vrs = lambda x: x.apply(lambda row: row['fecha_upc'] if row['VRS_D1'] == 1 else pd.NaT, axis=1),
              fecha_upc_vrs_Dall = lambda x: x.apply(lambda row: row['fecha_upc'] if row['VRS_Dall'] == 1 else pd.NaT, axis=1),
              
              
              #################################################################################################################################################
            #   fecha_upc_vrs = lambda x: x.fecha_upc_vrs_Dall,
              #################################################################################################################################################
              
              
              days_upc = 0,
              dias_en_ing = 0,
              days_estad_vrs = lambda x: np.where((x['VRS_D1']==1) & (x['DIAS_ESTAD']>0), x['DIAS_ESTAD'], 0),
              days_estad_vrs_Dall = lambda x: np.where((x['VRS_Dall']==1) & (x['DIAS_ESTAD']>0), x['DIAS_ESTAD'], 0),
              
              #Others covs
              is_rural = lambda x: np.where(x.porcent_rural>0.5,1,0),
              categori_macro = lambda x: x['Macrozona2'].astype('category').cat.codes + 1,
              categori_regions = lambda x: x['region'].astype('category').cat.codes + 1,
              exp_rural = lambda x: np.exp(x['porcent_rural']),
              percent_poor = lambda x: x['percent_poor'].fillna(x['percent_poor'].mean()) * 100,
              percent_poor_multidim = lambda x: x['percent_poor_multidim'].fillna(x['percent_poor_multidim'].mean()),
              is_poor = lambda x: np.where(x.percent_poor_multidim >= x.percent_poor_multidim.median(), 1, 0),
              vrs_pre_campaña_2024 = lambda x: np.where(x.fechaIng_vrs < pd.to_datetime("2024-04-01"), 1, 0), ################################### '2023-11-01'
              lrti_pre_campaña_2024  = lambda x: np.where(x.fechaIng_LRTI < pd.to_datetime("2024-04-01"), 1, 0),
              any_pre_campaña_2024  = lambda x: np.where(x.fechaIng_any < pd.to_datetime("2024-04-01"), 1, 0),
              upc_pre_campaña_2024  = lambda x: np.where(x.fecha_upc_vrs < pd.to_datetime("2024-04-01"), 1, 0),
              
              vrs_pre_campaña_2025 = lambda x: x.fechaIng_vrs.between(pd.to_datetime("2024-10-01"), pd.to_datetime("2025-03-01"), inclusive="neither").astype(int),
              lrti_pre_campaña_2025  = lambda x: x.fechaIng_vrs.between(pd.to_datetime("2024-10-01"), pd.to_datetime("2025-03-01"), inclusive="neither").astype(int),
              any_pre_campaña_2025  = lambda x: x.fechaIng_vrs.between(pd.to_datetime("2024-10-01"), pd.to_datetime("2025-03-01"), inclusive="neither").astype(int),
              upc_pre_campaña_2025  = lambda x: x.fechaIng_vrs.between(pd.to_datetime("2024-10-01"), pd.to_datetime("2025-03-01"), inclusive="neither").astype(int),
              )
    )

    print(np.unique(reg_print))
    
    # DAYS UPC
    for i in range(1,10):
        diff = f'dias_en_area_{i}'
        df[diff] = 0
    
    for index, row in df.copy().query('cama=="UPC"').iterrows():
        if pd.isna(row['fecha_tras_1']):
            row['days_upc'] = (row['fechaEgr'] - row['fechaIng_any']).days
            
        else:
            row['dias_en_ing'] = (row['fecha_tras_1'] - row['fechaIng_any']).days
        
            for i in range(1, 10):
                date_col = f'fecha_tras_{i}'
                date_col_next = f'fecha_tras_{i+1}'
                diff = f'dias_en_area_{i}'
                
                if date_col_next not in tras_date.values() or pd.isna(row[date_col_next]):
                    date_col_next = 'fechaEgr'
                    row[diff] = (row[date_col_next] - row[date_col]).days
                    break
                
                row[diff] = (row[date_col_next] - row[date_col]).days

            row['days_upc'] += row['dias_en_ing']*row['AREA_FUNC_I']
            
            for i in range(1, 10):
                area_col = f'AREAF_{i}_TRAS'
                diff = f'dias_en_area_{i}'
                if row[area_col]==1:
                    row['days_upc'] += row[diff]
            
        df.loc[index, 'days_upc'] = np.where(row['days_upc']==0, 1, row['days_upc'])
        df.loc[index, 'dias_en_ing'] = row['dias_en_ing']
        
        for i in range(1,10):
            diff = f'dias_en_area_{i}'
            df.loc[index, diff] = row[diff]

    df = df.assign(
        days_upc_vrs = lambda x: np.where((x['VRS_D1']==1) & (x['days_upc']>0), x['days_upc'], 0),
        days_upc_vrs_Dall = lambda x: np.where((x['VRS_Dall']==1) & (x['days_upc']>0), x['days_upc'], 0),
        
        #Events
        event_upc = lambda x : (x['fecha_upc_vrs'].notnull()),
        event_upc_Dall = lambda x : (x['fecha_upc_vrs_Dall'].notnull()),
        event_vrs = lambda x : x['fechaIng_vrs'].notnull(),
        event_vrs_Dall = lambda x : x['fechaIng_vrs_Dall'].notnull(),
        event_LRTI = lambda x: x['fechaIng_LRTI'].notnull(),
        event_any = lambda x : x['fechaIng_any'].notnull(),
        take_nirse = lambda x : x['fechaInm'].notnull()
        )
    print(f'n_rows_post_prefiltred= {df.shape[0]}')
    
    return df


def filtros_IH_new(df,
                   semana_lower=24,
                   semana_upper=42,
                   eliminar_inmunes_pre_season=True, 
                   T_inicial = pd.to_datetime('2025-03-01'), 
                   fecha_dt=None,
                   duracion_dias_nirse=180,
                   meses_inm=None,
                   ef_2024_in_2025=False,
                   fecha_cohort_in=None,
                   fecha_cohort_out=None):

    fecha_in = fecha_cohort_in
    fecha_out = fecha_cohort_out
    df = df.query('@fecha_in <= fecha_nac <= @fecha_out').copy()
    # df = df.query('(@fecha_in <= fechaInm <= @fecha_out) | (fechaInm.isna())').copy()  ## filtro para método 2 de duración

    ################################################################################################## ELIMINAR A INMUNIZADOS DSPS DE SEASON ###########################
    df = df[~(df.fechaInm >= fecha_dt)]  # es mejor usar negacion, pq si fechainm es na da true
    ###################################################################################################################################################################
    
    #MUERTOS
    if 'VIVO' in(df.columns):
        print("Datos perdidos por muertes: ", df.query('VIVO=="NO"').shape[0])
        df_vivos = df.copy().query('VIVO=="SI"')
        
    n_ruts_ini = df_vivos[(df_vivos['SEMANAS'].notna()) & (df_vivos['PESO'].notna()) & (df_vivos['p_00001_lognormal'].notna()) & (df_vivos['p_99999_lognormal'].notna())].RUN.nunique()
    
    print("ruts perdidos por filtro semanas y peso: " , df_vivos[((df_vivos['SEMANAS']<semana_lower) | (df_vivos['SEMANAS']>semana_upper) |
                                                                   (df_vivos['PESO'] < df_vivos['p_00001_lognormal']) | (df_vivos['PESO'] > df_vivos['p_99999_lognormal']))].RUN.nunique())
    
    #SEMANA Y PESO GESTACION
    df_vivos_f = df_vivos[(df_vivos['SEMANAS']>=semana_lower) & (df_vivos['SEMANAS']<=semana_upper) &
                      (df_vivos['PESO'] >= df_vivos['p_00001_lognormal']) & (df_vivos['PESO'] <= df_vivos['p_99999_lognormal'])]
    
   
    #drop intersex
    df1 = df_vivos_f.query('SEXO!=9')
    print('Droped intersex:', df_vivos_f.query('SEXO==9').shape[0])
    
    #EDAD MADRE ATIPICA
    print("Datos perdidos por edad madre atípica:", (df1.RUN.nunique() - df1[(df1['EDAD_M']<=51) & (df1['EDAD_M']>=12)].RUN.nunique()))
    df1 = df1[(df1['EDAD_M']<=51) & (df1['EDAD_M']>=12)]
    
    #FECHA ING < FECHA NAC
    df1_outIngNac = df1[(df1['fechaIng_any'] < df1['fecha_nac'])]
    print("Datos perdidos por fecha ingreso menor a fecha nacimiento:", df1[df1.RUN.isin(df1_outIngNac.RUN.unique())].RUN.nunique())
    df1 = df1[~df1.RUN.isin(df1_outIngNac.RUN.unique())]
    
    df1.loc[df1.fechaInm < pd.to_datetime("2024-04-01"), 'fechaInm'] = pd.to_datetime("2024-04-01")
        
    if eliminar_inmunes_pre_season:
        df1 = df1[~(df1.fechaInm.between(pd.to_datetime("2024-09-30"), pd.to_datetime("2025-03-01"), inclusive="neither"))]
    
    df1 = df1.assign(fechaIng_vrs_copy = lambda x: x.fechaIng_vrs)
    
    df1.loc[df1.fechaIng_any <= df1.fecha_nac + pd.DateOffset(days=7), ['fechaIng_any', 'fechaIng_vrs', 'fecha_upc_vrs','fechaIng_LRTI']] = pd.NaT
    
    # df1 = df1[~df1.RUN.isin(df1[df1.fechaIng_any < pd.to_datetime("2024-04-01")].RUN.unique())]
    
    if 'MARCA' in(df1.columns):
        df1 = df1.query('MARCA!=1') ###############################################################
    
    print('vrs en los primeros 7 dias de nacer:', df1[df1.RUN.isin(df1[df1['fechaIng_vrs_copy'] <= df1.fecha_nac + pd.DateOffset(days=7)].RUN.unique())].RUN.nunique())
    
    ###############
    df1 = df1[~df1.RUN.isin(df1[df1['fechaIng_vrs_copy'] <= df1.fecha_nac + pd.DateOffset(days=7)].RUN.unique())]
    ###############
    
    n_ruts_fin = df1.RUN.nunique()
    
    print('Ruts eliminados:', n_ruts_ini - n_ruts_fin)
    
    #FECHA INGRESO < FECHA INMUNE
    fechas_df = {'fechaIng_vrs': df1, 'fecha_upc_vrs': df1, 'fechaIng_LRTI': df1, 'fechaIng_any': df1}
    
    for key in fechas_df.keys():
        df_onedit = fechas_df[key].copy()
        if key == 'fecha_upc_vrs':
            col_fecha_ing = 'fechaIng_vrs'
        else:
            col_fecha_ing = key
        
        if ef_2024_in_2025:
            df_onedit = df_onedit[~df_onedit.RUN.isin(df_onedit[(df_onedit[col_fecha_ing] < T_inicial) &                 ############ EFF 2024 en 2025
                                                                (df_onedit[col_fecha_ing] >= pd.to_datetime("2024-07-01"))].RUN.unique())]
            
            if col_fecha_ing.endswith('vrs'):
                    df_onedit.loc[df_onedit[col_fecha_ing] < T_inicial, col_fecha_ing] = pd.NaT ############ EFF 2024 en 2025 
        else:
            # df_onedit = df_onedit[~df_onedit.RUN.isin(df_onedit[df_onedit['fechaIng_vrs'] < T_inicial].RUN.unique())]
            # df_onedit = df_onedit[~df_onedit.RUN.isin(df_onedit[(df_onedit['fechaIng_vrs'] < T_inicial) &
            #                                                     (df_onedit['fecha_nac'] > pd.to_datetime("2024-09-30"))].RUN.unique())]
            
            df_onedit = df_onedit[~df_onedit.RUN.isin(df_onedit[df_onedit[col_fecha_ing] < T_inicial].RUN.unique())]
            
        df_onedit.loc[(df_onedit['fechaInm'] > df_onedit[col_fecha_ing]) , 'fechaInm'] = pd.NaT 
        print(key, 'Reemplazos n/a net 7 days inmunizado: ', df_onedit.loc[((df_onedit[col_fecha_ing] - df_onedit['fechaInm']).dt.days <= 7)].RUN.nunique())
        
        ########################################
        df_onedit.loc[((df_onedit[col_fecha_ing] - df_onedit['fechaInm']).dt.days <= 7)  , 'fechaInm'] = pd.NaT  
        ########################################
        fechas_df[key] = df_onedit
        
    #BORRAR DUPLICADOS
    df_aux = fechas_df['fecha_upc_vrs'].copy()
    df_upc = df_aux.copy()[df_aux.fecha_upc_vrs.notna()].sort_values(by='fecha_upc_vrs').drop_duplicates(subset=['RUN'], keep='first')
    df_noupc = df_aux.copy()[~df_aux.RUN.isin(df_upc.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_upc = pd.concat([df_noupc,df_upc])
    
    df_aux = fechas_df['fechaIng_vrs'].copy()
    vrs = df_aux.copy()[df_aux.fechaIng_vrs.notna()].sort_values(by='fechaIng_vrs').drop_duplicates(subset=['RUN'], keep='first')
    df_novrs = df_aux.copy()[~df_aux.RUN.isin(vrs.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_vrs = pd.concat([df_novrs,vrs])
    
    df_aux = fechas_df['fechaIng_LRTI'].copy()
    LRTI = df_aux.copy()[df_aux.fechaIng_LRTI.notna()].sort_values(by='fechaIng_LRTI').drop_duplicates(subset=['RUN'], keep='first')
    df_noLRTI = df_aux.copy()[~df_aux.RUN.isin(LRTI.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_LRTI = pd.concat([df_noLRTI,LRTI])
    
    df_aux = fechas_df['fechaIng_any'].copy()
    df_any = df_aux.copy()[df_aux.fechaIng_any.notna()].sort_values(by='fechaIng_any').drop_duplicates(subset=['RUN'], keep='first')
    df_noany = df_aux.copy()[~df_aux.RUN.isin(df_any.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_any = pd.concat([df_noany,df_any])
    
    #Crear cumplemeses vivo y cumplesemana inmune
    
    fechas_events = {'fechaInm':'inmunizado', 'fecha_upc_vrs':'event_upc','fechaIng_vrs':'event_vrs',
                     'fechaIng_LRTI':'event_LRTI', 'fechaIng_any':'event_any'}
    
    dfs_fecha = {'fechaIng_vrs': df_filtrado_vrs, 'fecha_upc_vrs': df_filtrado_upc, 'fechaIng_LRTI': df_filtrado_LRTI, 'fechaIng_any': df_filtrado_any}
    
    dfs_fecha_actualiza = {}
    
    for key, df_f in dfs_fecha.items():
        
        for fecha_col, event_col in fechas_events.items():
            df_f.loc[df_f[fecha_col] >= fecha_dt, fecha_col] = pd.NaT           
            df_f[event_col] = df_f[fecha_col].notnull().astype(int)   

       # df_save = df_f.loc[(df_f.fecha_nac <= pd.to_datetime("2024-09-30")) | (df_f.fecha_nac.dt.isocalendar().year == 2023)] 
        df_save = df_f.copy()
        dfs_fecha_actualiza[key] = df_save
        
    df_filtrado_any = dfs_fecha_actualiza['fechaIng_any']
    df_filtrado_LRTI = dfs_fecha_actualiza['fechaIng_LRTI']
    df_filtrado_vrs = dfs_fecha_actualiza['fechaIng_vrs']
    df_filtrado_upc = dfs_fecha_actualiza['fecha_upc_vrs']
        
    return df_filtrado_any, df_filtrado_LRTI, df_filtrado_vrs, df_filtrado_upc

def filtros_IH_case_control(df,dias=dias, meses=meses, semana_lower=24,semana_upper=42):

    #MUERTOS
    # if 'VIVO' in(df.columns):
    #     print("Datos perdidos por muertes: ", df.query('VIVO=="NO"').shape[0])
    #     df_vivos = df.copy().query('VIVO=="SI"')
        
    df = df.assign(ingreso_mayor_1_mes = lambda x:  (x.fechaIng_any - x.fecha_nac).dt.days)
    
    df['VIVO_mix'] = np.where((df['VIVO']=='SI')  , 'SI', 'NO') #| (df['VIVO'].isna()) & (df['FECHA_DEF'].isna())
    
    df_vivos = df.query('(VIVO_mix=="SI") | ((VIVO_mix=="NO") & (ingreso_mayor_1_mes>=7))')
        
    n_ruts_ini = df_vivos[(df_vivos['SEMANAS'].notna()) & (df_vivos['PESO'].notna()) & (df_vivos['p_00001_lognormal'].notna()) & (df_vivos['p_99999_lognormal'].notna())].RUN.nunique()
    
    # print("ruts perdidos por filtro semanas y peso: " , df_vivos[((df_vivos['SEMANAS']<semana_lower) | (df_vivos['SEMANAS']>semana_upper) |
    #                                                                (df_vivos['PESO'] < df_vivos['p_00001_lognormal']) | (df_vivos['PESO'] > df_vivos['p_99999_lognormal']))].RUN.nunique())
    
    #SEMANA Y PESO GESTACION
    df_vivos_f = df_vivos[(df_vivos['SEMANAS']>=semana_lower) & (df_vivos['SEMANAS']<=semana_upper) &
                      (df_vivos['PESO'] >500) & (df_vivos['PESO'] <= df_vivos['p_99999_lognormal'])]
    
        # df_vivos_f = df_vivos[(df_vivos['SEMANAS']>=semana_lower) & (df_vivos['SEMANAS']<=semana_upper) &
        #               (df_vivos['PESO'] >= df_vivos['p_00001_lognormal']) & (df_vivos['PESO'] <= df_vivos['p_99999_lognormal'])]
    
    #df_vivos_f = df_vivos
    # print("Datos perdidos por filtro peso: " , df.drop_duplicates().shape[0]-df_filtered.drop_duplicates().shape[0])
    
    
    #drop intersex
    df1 = df_vivos_f.query('SEXO!=9')
    print('Droped intersex:', df_vivos_f.query('SEXO==9').shape[0])
    
    #EDAD MADRE ATIPICA
    print("Datos perdidos por edad madre atípica:", (df1.RUN.nunique() - df1[(df1['EDAD_M']<=51) & (df1['EDAD_M']>=12)].RUN.nunique()))
    df1 = df1[(df1['EDAD_M']<=51) & (df1['EDAD_M']>=12)]
    
    #FECHA ING < FECHA NAC
    df1_outIngNac = df1[(df1['fechaIng_any'] < df1['fecha_nac'])]
    print("Datos perdidos por fecha ingreso menor a fecha nacimiento:", df1[df1.RUN.isin(df1_outIngNac.RUN.unique())].RUN.nunique())
    df1 = df1[~df1.RUN.isin(df1_outIngNac.RUN.unique())]
    
    df1.loc[df1.fechaInm < pd.to_datetime("2024-04-01"), 'fechaInm'] = pd.to_datetime("2024-04-01")
    
    df1 = df1.assign(fechaIng_vrs_copy = lambda x: x.fechaIng_vrs)
    
    df1.loc[df1.fechaIng_any <= df1.fecha_nac + pd.DateOffset(days=7), ['fechaIng_any', 'fechaIng_vrs', 'fecha_upc_vrs','fechaIng_LRTI']] = pd.NaT
    
    # df1 = df1[~df1.RUN.isin(df1[df1.fechaIng_any < pd.to_datetime("2024-04-01")].RUN.unique())]
    
    if 'MARCA' in(df1.columns):
        df1 = df1.query('MARCA!=1') ###############################################################
    
    
    print('vrs en los primeros 7 dias de nacer:', df1[df1.RUN.isin(df1[df1['fechaIng_vrs_copy'] <= df1.fecha_nac + pd.DateOffset(days=7)].RUN.unique())].RUN.nunique())
    
    ###############
    df1 = df1[~df1.RUN.isin(df1[df1['fechaIng_vrs_copy'] <= df1.fecha_nac + pd.DateOffset(days=7)].RUN.unique())]
    ###############
    
    n_ruts_fin = df1.RUN.nunique()
    
    print('Ruts eliminados:', n_ruts_ini - n_ruts_fin)
    
    #FECHA INGRESO < FECHA INMUNE
    fechas_df = {'fechaIng_vrs': df1, 'fecha_upc_vrs': df1, 'fechaIng_LRTI': df1, 'fechaIng_any': df1}
    
    for key in fechas_df.keys():
        df_onedit = fechas_df[key].copy()
        if key == 'fecha_upc_vrs':
            df_onedit = df_onedit[~df_onedit.RUN.isin(df_onedit[df_onedit['fechaIng_vrs'] < pd.to_datetime("2024-04-01")].RUN.unique())]
        else:
            df_onedit = df_onedit[~df_onedit.RUN.isin(df_onedit[df_onedit[key] < pd.to_datetime("2024-04-01")].RUN.unique())]
        
        df_onedit.loc[(df_onedit['fechaInm'] > df_onedit[key]) , 'fechaInm'] = pd.NaT 
        print(key, 'Reemplazos n/a net 7 days inmunizado: ', df_onedit.loc[((df_onedit[key] - df_onedit['fechaInm']).dt.days <= 7)].RUN.nunique()) 
        ########################################
        df_onedit.loc[((df_onedit[key] - df_onedit['fechaInm']).dt.days <= 7)  , 'fechaInm'] = pd.NaT  
        ########################################
        fechas_df[key] = df_onedit
        
    #BORRAR DUPLICADOS
    df_aux = fechas_df['fecha_upc_vrs'].copy()
    df_upc = df_aux.copy()[df_aux.fecha_upc_vrs.notna()].sort_values(by='fecha_upc_vrs').drop_duplicates(subset=['RUN'], keep='first')
    df_noupc = df_aux.copy()[~df_aux.RUN.isin(df_upc.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_upc = pd.concat([df_noupc,df_upc])
    
    df_aux = fechas_df['fechaIng_vrs'].copy()
    vrs = df_aux.copy()[df_aux.fechaIng_vrs.notna()].sort_values(by='fechaIng_vrs').drop_duplicates(subset=['RUN'], keep='first')
    df_novrs = df_aux.copy()[~df_aux.RUN.isin(vrs.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_vrs = pd.concat([df_novrs,vrs])
    
    df_aux = fechas_df['fechaIng_LRTI'].copy()
    LRTI = df_aux.copy()[df_aux.fechaIng_LRTI.notna()].sort_values(by='fechaIng_LRTI').drop_duplicates(subset=['RUN'], keep='first')
    df_noLRTI = df_aux.copy()[~df_aux.RUN.isin(LRTI.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_LRTI = pd.concat([df_noLRTI,LRTI])
    
    df_aux = fechas_df['fechaIng_any'].copy()
    df_any = df_aux.copy()[df_aux.fechaIng_any.notna()].sort_values(by='fechaIng_any').drop_duplicates(subset=['RUN'], keep='first')
    df_noany = df_aux.copy()[~df_aux.RUN.isin(df_any.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_any = pd.concat([df_noany,df_any])
    
    #Crear cumplemeses vivo y cumplesemana inmune
    
    fechas_events = {'fechaInm':'inmunizado', 'fecha_upc_vrs':'event_upc','fechaIng_vrs':'event_vrs',
                     'fechaIng_LRTI':'event_LRTI', 'fechaIng_any':'event_any'}
    
    dfs_fecha = {'fechaIng_vrs': df_filtrado_vrs, 'fecha_upc_vrs': df_filtrado_upc, 'fechaIng_LRTI': df_filtrado_LRTI, 'fechaIng_any': df_filtrado_any}
    
    dfs_fecha_actualiza = {}
    
    fecha_dt = pd.to_datetime('2024-09-30') #+ pd.DateOffset(days=7) 2024-09-30
    
    for key, df_f in dfs_fecha.items():
    
        # for m in meses:
        #     nombre_columna = f'age_{m}m'
        #     nombre_columna_bin = f'si_{m}_meses'
        #     df_f[nombre_columna] = df_f['fecha_nac'] + pd.DateOffset(months=m)
        #     df_f.loc[(df_f[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
        #     df_f[nombre_columna_bin] = 1 - df_f[nombre_columna].isna()
        
        # for d in dias:
        #     nombre_columna = f'inm_{d}d'
        #     nombre_columna_bin = f'inm_mayor_{d}d'
        #     df_f[nombre_columna] = df_f['fechaInm'] + pd.DateOffset(days=d)
        #     df_f.loc[(df_f[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
        #     df_f[nombre_columna_bin] = 1 - df_f[nombre_columna].isna()
        
        for fecha_col, event_col in fechas_events.items():
            df_f.loc[df_f[fecha_col] > fecha_dt, fecha_col] = pd.NaT
            df_f[event_col] = df_f[fecha_col].notnull().astype(int)
    
        df_save = df_f.loc[(df_f.fecha_nac <= pd.to_datetime("2024-09-30")) | (df_f.fecha_nac.dt.isocalendar().year == 2023)] #2024-09-30

        dfs_fecha_actualiza[key] = df_save
        
    df_filtrado_any = dfs_fecha_actualiza['fechaIng_any']
    df_filtrado_LRTI = dfs_fecha_actualiza['fechaIng_LRTI']
    df_filtrado_vrs = dfs_fecha_actualiza['fechaIng_vrs']
    df_filtrado_upc = dfs_fecha_actualiza['fecha_upc_vrs']
        
    return df_filtrado_any, df_filtrado_LRTI, df_filtrado_vrs, df_filtrado_upc

def call_data(path=None):
    
    leo_zonas_rename = {'Macrozona Centro Sur':'South Macrozone',
                        'Macrozona Sur':'South Macrozone',
                        'Macrozona Norte':'North Macrozone',
                        'Macrozona Centro Norte':'Central Macrozone',
                        'Macrozona Austral':'Austral Macrozone'}
    region_to_macrozone_agencia = {
        "ARICA Y PARINACOTA": "Macrozona Norte",
        "TARAPACA": "Macrozona Norte",
        "ANTOFAGASTA": "Macrozona Norte",
        "ATACAMA": "Macrozona Norte",
        "COQUIMBO": "Macrozona Centro Norte",
        "VALPARAISO": "Macrozona Centro Norte",
        "METROPOLITANA": "Macrozona Centro Norte",
        "O'HIGGINS": "Macrozona Centro Norte",
        "MAULE": "Macrozona Centro Sur",
        "NUBLE": "Macrozona Centro Sur",
        "BIOBIO": "Macrozona Centro Sur",
        "ARAUCANIA": "Macrozona Centro Sur",
        "LOS RIOS": "Macrozona Sur",
        "LOS LAGOS": "Macrozona Sur",
        "AISEN": "Macrozona Austral",
        "MAGALLANES Y ANTARTICA": "Macrozona Austral"
    }

    path_actual = Path.cwd()
    path_data = path_actual.parent / 'Data'

    #ruts_cariopaticos_1 = pd.read_csv(path_data/'ruts_cardiopatias.csv', index_col=0).query('card1')
    #ruts_cariopaticos_2 = pd.read_csv(path_data/'ruts_cardiopatias.csv', index_col=0).query('card2')
    
    with open(path_data/'lista_ruts_cardio.pkl', 'rb') as f:
        lista_ruts_cardio = pickle.load(f)
    
    df_pf = (
        pre_filtred(df_name=path)
        # Drop de columnas según nombre
        .pipe(lambda df: df.drop(columns=[
            col for col in df.columns
            if (
                col.startswith('inm_') or
                col.startswith('MES_') or
                col.startswith('DIA_') or
                col.startswith('SERC_') or
                (col.startswith('ANO_') and col.endswith('TRAS')) or
                col.startswith('AREAF') or
                col.startswith('dias_en') or
                col.startswith('fecha_tras') or
                col.endswith('Dall')
            )
        ]))

        .pipe(lambda df: df.assign(
            region=df.region.where(df.RUN != 'ac483764636448d753930868dd3192a785e3728a2d81b3fc44dba45c3506255a', 'METROPOLITANA')
        ))
        .pipe(lambda df: df.assign(
            region=df.region.where(df.COMUNA != 12202, 'MAGALLANES Y ANTARTICA')
        ))

        .assign(
            cardio1=lambda x: x.RUN.isin(lista_ruts_cardio).astype(int),
            log_peso=lambda x: np.log(x.PESO),
            Macrozones = lambda x: x.region.map(region_to_macrozone_agencia).replace(leo_zonas_rename),
            super_preterm = lambda x: ((x.SEMANAS<32) | (x.PESO<1500)).astype(int),
        )
        .copy()
    )

    df_f_any, df_f_LRTI, df_f_vrs, df_f_upc = filtros_IH_case_control(df_pf)

    df_f_vrs_cut = df_f_vrs.copy()
    df_f_upc_cut = df_f_upc.copy()
    
    df_upc = df_f_upc_cut.copy().rename(columns={'event_upc':'event'})
    df_vrs = df_f_vrs_cut.copy().rename(columns={'event_vrs':'event'})
    
    df_dic = {'vrs':df_vrs, 'upc':df_upc}
    df_return = {}
    
    for key, df in df_dic.items():
        
        # df_case = (df
        #         .query("event == True")
        #         #.drop_duplicates(subset=['RUN'], keep='first') #innecesario, no hay duplicados
        #         .copy()
        #         )

        # # Paso 3: De ese subset, quedarnos con el primer caso de cada 'RUN' donde 'diagnostico1' no es nulo
        # df_non_vrs = (df[~df['RUN'].isin(df_case['RUN'])]
        #             #.drop_duplicates(subset=['RUN'], keep='first')
        #             .query("DIAG1.notnull()")
        #             )

        # # Paso 4: Filtrar los casos donde 'diagnostico1' es nulo en el subset original sin VRS
        # df_nulo_diag = (df[~df['RUN'].isin(df_case['RUN']) & ~df['RUN'].isin(df_non_vrs['RUN'])]
        #             .query("DIAG1.isnull()")
        #             #.drop_duplicates(subset=['RUN'], keep='first')
        #             )

        fecha_referencia_nacido_str = '2024-09-30'
        fecha_referencia_campaña_str = '2024-09-30'
        chile_chico = ['METROPOLITANA','VALPARAISO','LOS LAGOS','LOS RIOS','MAULE','TARAPACA']

        df_final = (#pd.concat([df_case, df_non_vrs, df_nulo_diag], ignore_index=True)
                    df
                    .assign(
                    chile_chico = lambda df: np.where(df['region'].isin(chile_chico),1,0),
                    estado_inmunizacion = lambda df: df['inmunizado'],#.map({'inmunizado': 1, 'not inmunizado': 0})
                    diag_vrs =  lambda df: df['event'].astype(int),
                    diag_1_leter = lambda x: np.where(x.fechaIng_any.notna(), x['DIAG1'].str[0], pd.NA),
                    edad_relativa =  lambda x: (x['fecha_nac']-pd.to_datetime(fecha_referencia_nacido_str)).dt.days,
                    ingreso_relativo =  lambda x: (x['fechaIng_any']-pd.to_datetime(fecha_referencia_campaña_str)).dt.days
                    )
                    
        )
        
        df_return[key] = df_final.copy()
        
    return df_return['vrs'], df_return['upc']

def call_data_2025_update(path=None,
                          fecha_referencia_nacido_str = '2024-10-01',
                          fecha_referencia_campaña_str = '2025-03-01',
                          fecha_cohort_in = pd.to_datetime('2024-10-01'),
                          fecha_cohort_out = pd.to_datetime('2025-09-30'),
                          T_inicial = pd.to_datetime('2025-03-01'), 
                          fecha_dt = None,
                          eliminar_inmunes_pre_season = True,
                          ef_2024_in_2025 = True):
    
    leo_zonas_rename = {'Macrozona Centro Sur':'South Macrozone',
                        'Macrozona Sur':'South Macrozone',
                        'Macrozona Norte':'North Macrozone',
                        'Macrozona Centro Norte':'Central Macrozone',
                        'Macrozona Austral':'Austral Macrozone'}
    
    region_to_macrozone_agencia = {
        "ARICA Y PARINACOTA": "Macrozona Norte",
        "TARAPACA": "Macrozona Norte",
        "ANTOFAGASTA": "Macrozona Norte",
        "ATACAMA": "Macrozona Norte",
        "COQUIMBO": "Macrozona Centro Norte",
        "VALPARAISO": "Macrozona Centro Norte",
        "METROPOLITANA": "Macrozona Centro Norte",
        "O'HIGGINS": "Macrozona Centro Norte",
        "MAULE": "Macrozona Centro Sur",
        "NUBLE": "Macrozona Centro Sur",
        "BIOBIO": "Macrozona Centro Sur",
        "ARAUCANIA": "Macrozona Centro Sur",
        "LOS RIOS": "Macrozona Sur",
        "LOS LAGOS": "Macrozona Sur",
        "AISEN": "Macrozona Austral",
        "MAGALLANES Y ANTARTICA": "Macrozona Austral"
    }

    path_actual = Path.cwd()
    path_data = path_actual.parent / 'Data'
    
    with open(path_data/'lista_ruts_cardio.pkl', 'rb') as f:
        lista_ruts_cardio = pickle.load(f)
    
    df_pf = (
        pre_filtred(df_name=path)
        # Drop de columnas según nombre
        .pipe(lambda df: df.drop(columns=[
            col for col in df.columns
            if (
                col.startswith('inm_') or
                col.startswith('MES_') or
                col.startswith('DIA_') or
                col.startswith('SERC_') or
                (col.startswith('ANO_') and col.endswith('TRAS')) or
                col.startswith('AREAF') or
                col.startswith('dias_en') or
                col.startswith('fecha_tras') or
                col.endswith('Dall')
            )
        ]))

        .pipe(lambda df: df.assign(
            region=df.region.where(df.RUN != 'ac483764636448d753930868dd3192a785e3728a2d81b3fc44dba45c3506255a', 'METROPOLITANA')
        ))
        .pipe(lambda df: df.assign(
            region=df.region.where(df.COMUNA != 12202, 'MAGALLANES Y ANTARTICA')
        ))

        .assign(
            cardio1=lambda x: x.RUN.isin(lista_ruts_cardio).astype(int),
            Macrozones = lambda x: x.region.map(region_to_macrozone_agencia).replace(leo_zonas_rename),
            super_preterm = lambda x: ((x.SEMANAS<32) | (x.PESO<1500)).astype(int),
        )
        .copy()
    )
    _, _, df_f_vrs, df_f_upc = filtros_IH_new(df_pf,
                                              semana_lower=24,
                                              semana_upper=42,
                                              eliminar_inmunes_pre_season=eliminar_inmunes_pre_season, 
                                              T_inicial = T_inicial, 
                                              fecha_dt=fecha_dt,
                                              ef_2024_in_2025=ef_2024_in_2025,
                                              fecha_cohort_in=fecha_cohort_in,
                                              fecha_cohort_out=fecha_cohort_out)
    df_f_vrs_cut = df_f_vrs.copy()
    df_f_upc_cut = df_f_upc.copy()
    
    df_upc = df_f_upc_cut.copy().rename(columns={'event_upc':'event'})
    df_vrs = df_f_vrs_cut.copy().rename(columns={'event_vrs':'event'})
    
    df_dic = {'vrs':df_vrs, 'upc':df_upc}
    df_return = {}
    
    for key, df in df_dic.items():
        
        fecha_referencia_nacido_str = '2023-10-28'
        fecha_referencia_campaña_str = '2024-04-01'
        chile_chico = ['METROPOLITANA','VALPARAISO','LOS LAGOS','LOS RIOS','MAULE','TARAPACA']

        df_final = (#pd.concat([df_case, df_non_vrs, df_nulo_diag], ignore_index=True)
                    df
                    .assign(
                    chile_chico = lambda df: np.where(df['region'].isin(chile_chico),1,0),
                    estado_inmunizacion = lambda df: df['inmunizado'],#.map({'inmunizado': 1, 'not inmunizado': 0})
                    diag_vrs =  lambda df: df['event'].astype(int),
                    diag_1_leter = lambda x: np.where(x.fechaIng_any.notna(), x['DIAG1'].str[0], pd.NA),
                    edad_relativa =  lambda x: (x['fecha_nac']-pd.to_datetime(fecha_referencia_nacido_str)).dt.days,
                    ingreso_relativo =  lambda x: (x['fechaIng_any']-pd.to_datetime(fecha_referencia_campaña_str)).dt.days
                    )

        )
        
        df_return[key] = df_final.copy()
        
    return df_return['vrs'], df_return['upc']

def flexible_matching(df_cases,
                      df_control,
                      match_vars,
                      match_date_vars=None,
                      intervals=None,
                      n_control=1
                      ):
    cases = df_cases.copy().reset_index(drop=True)
    controls = df_control.copy().reset_index(drop=True)
    
    matched_pairs = []
    run_attached= []
    cases_matched =[]
    contador=0
    
    # Para cada caso, buscamos controles que coincidan en las variables de matching
    for i, case in cases.iterrows():
        matched_controls = controls[~controls.RUN.isin(run_attached)].copy()  # Copiamos los controles disponibles

        # Recorremos las variables de emparejamiento y aplicamos reglas
        for var in match_vars:
            # Comprobamos si intervals no es None y contiene la variable
            if intervals and var in intervals:
                # Si hay un intervalo definido para esta variable, aplicamos el filtro
                matched_controls = matched_controls[
                    (matched_controls[var] >= (case[var] - intervals[var])) &
                    (matched_controls[var] <= (case[var] + intervals[var]))
                ]
            else:
                # Si no hay intervalo, hacemos un match exacto
                matched_controls = matched_controls[matched_controls[var] == case[var]]

        # Si se define una variable temporal para comparar fechas, aplicamos el filtro de ±15 días
        if match_date_vars:
            for date_var, date_intervals in match_date_vars.items():
                matched_controls = matched_controls[
                    (matched_controls[date_var] >= case[date_var] - pd.Timedelta(days=date_intervals)) &
                    (matched_controls[date_var] <= case[date_var] + pd.Timedelta(days=date_intervals))
                ]
    
        # Si encontramos controles emparejados, los agregamos a la lista de pares
        if not matched_controls.empty:
            #contador+=1
            #(contador,i+1)
            # Tomamos el primer control disponible para este caso
            selected_control = matched_controls.sample(n=min(n_control, len(matched_controls)), random_state=42)
            selected_control['Matched_Case_RUN'] = case['RUN']

            # Añadimos el control emparejado a la lista de pares
            matched_pairs.append(selected_control)

            run_attached+=list(selected_control.RUN.unique())
            cases_matched.append(case['RUN'])

    print(f'Total cases matched is : {len(cases_matched)}')
    # Concatenamos todos los controles emparejados
    matched_controls = pd.concat(matched_pairs, axis=0)
    matched_data = pd.concat([cases[cases.RUN.isin(cases_matched)], matched_controls], axis=0)
    matched_data['Group'] = matched_data['Matched_Case_RUN'].fillna(matched_data['RUN'])


    # Devolvemos los casos con sus controles emparejados
    #matched_data = pd.concat([cases, matched_controls], axis=0)
    return matched_data,matched_controls

def integer_programming_matching_gurobi(df_cases, df_control, match_vars, intervals, match_date_vars, n_control=1):
    
    cases = df_cases.copy().reset_index(drop=True)
    controls = df_control.copy().reset_index(drop=True)
    
    model = Model("Maximize_Matching")
    model.setParam('OutputFlag', 0)  
    
    feasible_controls = {}
    star_time = time.time()
    for i, case in cases.iterrows():

        controls_copy = df_control.copy().reset_index(drop=True)

        for var in match_vars:
            if intervals and var in intervals:
                controls_copy = controls_copy[
                    (controls_copy[var] >= (case[var] - intervals[var])) & 
                    (controls_copy[var] <= (case[var] + intervals[var]))
                ]
            else:
                controls_copy = controls_copy[controls_copy[var] == case[var]]
        
        if match_date_vars:
            for date_var, date_interval in match_date_vars.items():
                controls_copy = controls_copy[
                    (controls_copy[date_var] >= case[date_var] - pd.Timedelta(days=date_interval)) &
                    (controls_copy[date_var] <= case[date_var] + pd.Timedelta(days=date_interval))
                ]

        # Solo añadir controles que cumplen con las condiciones de matching
        feasible_controls[i] = controls_copy.index.tolist() if not controls_copy.empty else []
    end_time = time.time()
    print(f"creacion conjuntos A_i time: {end_time - star_time}")
    
    star_time = time.time()
    x = {} 
    for i, control_indices in feasible_controls.items():
        for j in control_indices:
            x[(i, j)] = model.addVar(vtype=GRB.BINARY, name=f"x_{i}_{j}")   
    end_time = time.time()
    print(f"creacion variables time: {end_time - star_time}")
    
    model.setObjective(sum(x[(i, j)] for (i, j) in x.keys()), GRB.MAXIMIZE)

    for i in feasible_controls:
        model.addConstr(sum(x[(i, j)] for j in feasible_controls[i]) <= n_control, name=f"max_controls_case_{i}")

    for j in controls_copy.index:
        model.addConstr(sum(x[(i, j)] for i in feasible_controls if j in feasible_controls[i]) <= 1, name=f"max_case_control_{j}")


    star_time = time.time()
    model.optimize()
    
    end_time = time.time()
    print(f"optimize model time: {end_time - star_time}")
    
    cases_m =[]
    matched_pairs = []
    
    star_time = time.time()
    for (i, j), var in x.items():
        if var.X > 0.5:  # Solo emparejar si la variable de decisión es 1
                selected_control = controls.iloc[[j]].copy()
                selected_control['Matched_Case_RUN'] = cases.loc[i, 'RUN']
                matched_pairs.append(selected_control)
                cases_m.append(cases.loc[i, 'RUN'])
                
    end_time = time.time()
    print(f"matched_data time: {end_time - star_time}")    
            
    # Construir el DataFrame `matched_data` en el mismo formato
    matched_controls = pd.concat(matched_pairs, axis=0) if matched_pairs else pd.DataFrame()
    cases_matched = cases[cases['RUN'].isin(matched_controls['Matched_Case_RUN'])]
    matched_data = pd.concat([cases_matched, matched_controls], axis=0)
    matched_data['Group'] = matched_data['Matched_Case_RUN'].fillna(matched_data['RUN'])
    
    print(f'Total cases matched is : {len(cases_m)}')
    
    return matched_data, feasible_controls

def match_nn_max_dist_weigths(df_control, 
                              df_case, 
                              match_vars_nn, 
                              match_vars_exact=[], 
                              match_vars_onehot=[], 
                              ratio="1:1", 
                              max_distance=10000, 
                              weights={},
                              with_replacement=False,
                              random_state=123):
    print(ratio, ratio.split(":"))
    rng = np.random.default_rng(random_state)
    controls = df_control.copy().reset_index(drop=True)
    cases = df_case.copy().reset_index(drop=True)
    
    def apply_weights(df, columns, weights):
        weighted_df = df.copy()
        for col in columns:
            weighted_df[col] *= weights.get(col, 1)  # Aplica el peso; si no hay peso, usa 1
        return weighted_df
        
    for var in match_vars_onehot:
        if pd.api.types.is_numeric_dtype(controls[var]):
            controls[var] = controls[var].astype(str)
            cases[var] = cases[var].astype(str)

    case_count, control_count = map(int, ratio.split(":"))
    matched_pairs = []
    matched_incompleto = []
    used_indices = set()  

    def filter_exact(df, row, match_vars_exact):
        for var in match_vars_exact:
            df = df[df[var] == row[var]]
        return df

    covariables_case_copy = cases[match_vars_nn].copy()
    covariables_control_copy = controls[match_vars_nn].copy()
    
        # Concatenar los casos y controles para ajustar el escalador en ambos conjuntos
    combined_data = pd.concat([covariables_case_copy, covariables_control_copy])

    # Ajustar el escalador en el conjunto combinado
    scaler = StandardScaler()
    scaler.fit(combined_data)

    covariables_case_scaled = pd.DataFrame(scaler.transform(covariables_case_copy), columns=match_vars_nn)
    covariables_control_scaled = pd.DataFrame(scaler.transform(covariables_control_copy), columns=match_vars_nn)

    covariables_case = apply_weights(covariables_case_scaled, match_vars_nn, weights)
    covariables_control = apply_weights(covariables_control_scaled, match_vars_nn, weights)
        

    for col in match_vars_onehot:
        encoder = OneHotEncoder(sparse=False)
        combined_col = pd.concat([cases[[col]], controls[[col]]]).drop_duplicates().reset_index(drop=True)
        encoder.fit(combined_col)
        
        encoded_col_case = encoder.transform(cases[[col]])
        encoded_col_control = encoder.transform(controls[[col]])
        
        onehot_columns = encoder.get_feature_names_out([col])
        encoded_col_case_df = pd.DataFrame(encoded_col_case, columns=onehot_columns)
        encoded_col_control_df = pd.DataFrame(encoded_col_control, columns=onehot_columns)
        
        covariables_case = pd.concat([covariables_case, encoded_col_case_df], axis=1)
        covariables_control = pd.concat([covariables_control, encoded_col_control_df], axis=1)

    
    if  int(ratio.split(":")[0]) <= int(ratio.split(":")[1]): #ratio == "1:1" or ratio == "1:2" or ratio == "1:3" or
        neighbors_count = control_count

        for i, case_row in cases.iterrows():
            potential_controls = filter_exact(controls, case_row, match_vars_exact)
            if potential_controls.empty:
                matched_incompleto.append(case_row['RUN'])
                continue

            covariables_control_filtered = covariables_control.loc[potential_controls.index]

            nn = NearestNeighbors(n_neighbors=len(covariables_control_filtered), algorithm='auto')
            nn.fit(covariables_control_filtered)

            distances, indices = nn.kneighbors([covariables_case.iloc[i]])

            # selected_controls = []
            # for idx, distance in zip(indices[0], distances[0]):
            #     control_idx = potential_controls.index[idx]

            #     # --- regla de disponibilidad según reemplazo ---
            #     available = True if with_replacement else (control_idx not in used_indices)

            #     if distance <= max_distance and available:
            #         selected_controls.append(controls.loc[[control_idx]].copy())

            #         # solo marcamos como "usado" si es sin reemplazo
            #         if not with_replacement:
            #             used_indices.add(control_idx)

            #     if len(selected_controls) == neighbors_count:
            #         break

            # if len(selected_controls) == neighbors_count:
            #     for control in selected_controls:
            #         control['Matched_Case_RUN'] = case_row['RUN']
            #         matched_pairs.append(control)
            # else:
            #     matched_incompleto.append(case_row['RUN'])
            
            
            # selected_controls = []

            dist_arr = distances[0]
            idx_arr = indices[0]

            candidate_info = []
            for idx, distance in zip(idx_arr, dist_arr):
                control_idx = potential_controls.index[idx]
                available = True if with_replacement else (control_idx not in used_indices)

                if distance <= max_distance and available:
                    candidate_info.append((control_idx, distance))

            if len(candidate_info) < neighbors_count:
                matched_incompleto.append(case_row['RUN'])
                continue

            # ordenar por distancia
            candidate_info.sort(key=lambda x: x[1])

            chosen_indices = []
            pos = 0

            while len(chosen_indices) < neighbors_count and pos < len(candidate_info):
                current_dist = candidate_info[pos][1]

                tied_group = []
                while pos < len(candidate_info) and np.isclose(candidate_info[pos][1], current_dist):
                    tied_group.append(candidate_info[pos][0])
                    pos += 1

                remaining = neighbors_count - len(chosen_indices)

                if len(tied_group) <= remaining:
                    chosen_indices.extend(tied_group)
                else:
                    sampled = rng.choice(tied_group, size=remaining, replace=False)
                    chosen_indices.extend(sampled.tolist())

            if len(chosen_indices) == neighbors_count:
                for control_idx in chosen_indices:
                    control = controls.loc[[control_idx]].copy()
                    control['Matched_Case_RUN'] = case_row['RUN']
                    matched_pairs.append(control)

                if not with_replacement:
                    used_indices.update(chosen_indices)
            else:
                matched_incompleto.append(case_row['RUN'])


        # for i, case_row in cases.iterrows():
        #     potential_controls = filter_exact(controls, case_row, match_vars_exact)
        #     if potential_controls.empty:
        #         matched_incompleto.append(case_row['RUN'])
        #         continue
                
        #     covariables_control_filtered = covariables_control.loc[potential_controls.index]

        #     # nn = NearestNeighbors(n_neighbors=min(neighbors_count*5, len(covariables_control_filtered)), algorithm='auto')
        #     nn = NearestNeighbors(n_neighbors=len(covariables_control_filtered), algorithm='auto')

        #     nn.fit(covariables_control_filtered)

        #     distances, indices = nn.kneighbors([covariables_case.iloc[i]])
        #     selected_controls = []
        #     for idx, distance in zip(indices[0], distances[0]):
        #         control_idx = potential_controls.index[idx]
        #         if distance <= max_distance and control_idx not in used_indices:
        #             selected_controls.append(controls.loc[[control_idx]].copy())
        #             used_indices.add(control_idx)
        #         if len(selected_controls) == neighbors_count:
        #             break

        #     if len(selected_controls) == neighbors_count:
        #         for control in selected_controls:
        #             control['Matched_Case_RUN'] = case_row['RUN']
        #             matched_pairs.append(control)
        #     else:
        #         matched_incompleto.append(case_row['RUN'])
        
        if matched_pairs==[]:
            print('here?????????????')
            return "No hay controles disponibles para emparejar"
        else:
            print('here aaa')
            matched_controls = pd.concat(matched_pairs, axis=0) if matched_pairs else pd.DataFrame()
            cases_matched = cases[cases['RUN'].isin(matched_controls['Matched_Case_RUN'])]
            matched_data = pd.concat([cases_matched, matched_controls], axis=0).reset_index(drop=True)
            matched_data['Group'] = matched_data['Matched_Case_RUN'].fillna(matched_data['RUN'])
            cases_m = cases_matched['RUN'].tolist()
            
            print(f'Total cases = {cases.shape[0]}, Total controls = {controls.shape[0]}')
            print(f'Total cases matched is : {len(cases_m)}, Total control matched is : {matched_controls.shape[0]}') 
            print('ratio: ' + ratio) 
            print(f'No matched : {len(matched_incompleto)}')
            
            n_case_matched = len(cases_m)
            n_control_matched = matched_controls.shape[0]

    else: #ratio == "2:1" or ratio == "3:1"
        neighbors_count = case_count  

        for j, control_row in controls.iterrows():
            potential_cases = filter_exact(cases, control_row, match_vars_exact)
            if potential_cases.empty:
                matched_incompleto.append(control_row['RUN'])
                continue
            
            covariables_case_filtered = covariables_case.loc[potential_cases.index]

            nn = NearestNeighbors(n_neighbors=min(neighbors_count*5, len(covariables_case_filtered)), algorithm='auto')
            nn.fit(covariables_case_filtered)

            distances, indices = nn.kneighbors([covariables_control.iloc[j]])

            selected_cases = []
            for idx, distance in zip(indices[0], distances[0]):
                case_idx = potential_cases.index[idx]
                if distance <= max_distance and case_idx not in used_indices:
                    selected_cases.append(cases.loc[[case_idx]].copy())
                    used_indices.add(case_idx)
                if len(selected_cases) == neighbors_count:
                    break

            if len(selected_cases) == neighbors_count:
                for case in selected_cases:
                    case['Matched_Case_RUN'] = control_row['RUN']
                    matched_pairs.append(case)
            else:
                matched_incompleto.append(control_row['RUN'])
        print('here ofcourse')
        matched_controls = pd.concat(matched_pairs, axis=0) if matched_pairs else pd.DataFrame()
        cases_matched = controls[controls['RUN'].isin(matched_controls['Matched_Case_RUN'])]
        matched_data = pd.concat([cases_matched, matched_controls], axis=0).reset_index(drop=True)
        matched_data['Group'] = matched_data['Matched_Case_RUN'].fillna(matched_data['RUN'])
        cases_m = cases_matched['RUN'].tolist()
        
        print(f'Total cases = {cases.shape[0]}, Total controls = {controls.shape[0]}')
        print(f'Total cases matched is : {matched_controls.shape[0]}, Total control matched is : {len(cases_m)}') 
        print('ratio: ' + ratio) 
        print(f'No matched : {len(matched_incompleto)}')
        
        n_control_matched = len(cases_m)
        n_case_matched = matched_controls.shape[0]
    
    return matched_data, matched_incompleto, matched_controls, cases, controls, ratio, n_control_matched, n_case_matched


def summary_inicial(df):
    #cuantos diagnosticos x letra
    print("cuantos diagnosticos x letra tengo?\n")
    print(f"Total Diag: {df[df.DIAG1.notnull()].shape[0]}")
    print((df
    .assign(diag_1_leter=lambda x: x['DIAG1'].str[0])
    .groupby(['diag_1_leter']).size().sort_values(ascending=False)
    .head(10)
    ))
    contingency_table = pd.crosstab(
            df['diag_vrs'], 
            df['estado_inmunizacion'], 
            rownames=['Diagnóstico VRS'], 
            colnames=['Estado de Inmunización'], 
            margins=True  # Incluye totales
    )
    # Mostrar la tabla de contingencia
    percentage_table = contingency_table.div(contingency_table.loc['All'], axis=1) * 100
    percentage_by_column = contingency_table.div(contingency_table.sum(axis=0), axis=1) * 100

    # Mostrar la tabla de contingencia
    print("Tabla de Contingencia:\n", contingency_table)

    # Mostrar la tabla de porcentajes
    print("\nTabla de Porcentajes:\n", percentage_table)
    #print("\nTabla de Porcentajes por Columna (%):\n", percentage_by_column)
    
    # Contar RUN únicos en total
    total_run_unicos = df['RUN'].nunique()

    # Contar RUN únicos donde estado_inmunizacion es 1
    run_unicos_inmunizados = df[df['estado_inmunizacion'] == 1]['RUN'].nunique()

    print(f'Total RUN únicos: {total_run_unicos}')
    print(f'RUN únicos con estado_inmunizacion=1: {run_unicos_inmunizados}')
    print(f'% RUN únicos con estado_inmunizacion=1: {run_unicos_inmunizados/total_run_unicos:.1%}')
    
    # Suponiendo que tu DataFrame se llama 'df'

    # 1. Cantidad total de RUN únicos por ESTAB
    total_run_por_estab = df.groupby('ESTAB')['RUN'].nunique().reset_index(name='total_RUN_unicos')

    # 2. Cantidad de RUN únicos con estado_inmunizacion == 1 por ESTAB
    run_inmunizados_por_estab = (
        df[df['estado_inmunizacion'] == 1]
        .groupby('ESTAB')['RUN']
        .nunique()
        .reset_index(name='RUN_inmunizados_unicos')
    )

    # 3. Combinar ambos resultados
    inumization_by_estab = (total_run_por_estab
              .merge(run_inmunizados_por_estab, on='ESTAB', how='left')
              .assign(RUN_inmunizados_unicos = lambda df: df['RUN_inmunizados_unicos'].fillna(0).astype(int))
              .assign(porc_RUN_inmunizados_unicos =lambda df: df['RUN_inmunizados_unicos']/df['total_RUN_unicos'])
              .sort_values(by='porc_RUN_inmunizados_unicos', ascending=False)
    )

    print(inumization_by_estab)
    
    df_case= df.query("diag_vrs==True").copy()
    df_control= df.query("diag_vrs==False").copy()
    list_ESTAB_case = df_case.ESTAB.unique()
    case_grouped_by_estab = df_case.groupby(['ESTAB', 'estado_inmunizacion'])

    
    # Probabilidad de estar inmunizado por establecimiento en los casos que no son vrs
    prob_inmunizado_por_estab_control = (df_control[df_control['ESTAB'].isin(list_ESTAB_case)]
                                 .groupby('ESTAB')
                                 .apply( lambda x: x['estado_inmunizacion'].mean())
    )
    prob_inmunizado_por_estab_case_inmunization = (df_control[df_control['ESTAB'].isin(list_ESTAB_case)]
                                 .groupby('ESTAB')
                                 .apply( lambda x: x['estado_inmunizacion'].mean())
    )
    
    # Agrupar los casos por ESTAB y estado_inmunizacion, y calcular el promedio de inmunización en controles
    result_df = (
        df_case.groupby(['ESTAB', 'estado_inmunizacion'])
        .apply(lambda group: df_control[df_control['ESTAB'] == group.name[0]]['estado_inmunizacion'].mean())
        .reset_index(name='inmunizacion_promedio_control')
    )

    # Mostrar el resultado
    print(f"Que tan probable es escoger alguien inmunizado en la base case: {df_case['estado_inmunizacion'].mean().round(2):.0%}")
    print(f"Que tan probable es escoger alguien inmunizado en la base control dado que tu match era inmunizado o no?")

    print(result_df.groupby(['estado_inmunizacion'])['inmunizacion_promedio_control'].mean())

def comparar_medias_test(df1, df2, columnas):
    
    # Función para convertir columnas específicas a variables dummy
    def convertir_a_dummy(df, columnas):
        # Convertir solo las columnas especificadas a variables dummy
        cols_dummies = df[columnas].select_dtypes(include=['object']).columns

        df_dummies = pd.get_dummies(df[columnas], columns=cols_dummies, drop_first=True)
        return df_dummies
    
    # Convertir las columnas específicas a variables dummy
    df1 = convertir_a_dummy(df1, columnas)
    df2 = convertir_a_dummy(df2, columnas)
    
    
    # Actualizar la lista de columnas tras la conversión
    columnas_actualizadas = [col for col in df1.columns if col.startswith(tuple(columnas))]
    #columnas_actualizadas = [col for col in columnas if col in df1.columns] + nuevas_columnas_dummies
    
    resultados = []
    
    for col in columnas_actualizadas:
        # Calcular medias de las columnas seleccionadas en cada DataFrame
        media_df1 = df1[col].mean()
        media_df2 = df2[col].mean()
        
        # Aplicar t-test de dos muestras
        t_stat, p_value = stats.ttest_ind(df1[col], df2[col], nan_policy='omit')
        # Añadir el mensaje si el p-valor es menor a 0.05
        mensaje = ""
        if p_value < 0.05:
            mensaje = "Existe una diferencia significativa."
        

        # Almacenar resultados en una lista de diccionarios
        resultados.append({
            'Columna': col,
            'Media_df1': media_df1,
            'Media_df2': media_df2,
            'T-stat': t_stat,
            'P-value': p_value,
            'Mensaje': mensaje  # Añadir el mensaje al resultado

        })
    
    # Convertir los resultados en un DataFrame para mejor visualización
    resultados_df = pd.DataFrame(resultados).round(2)
    return resultados_df

def analyze_vrs_data(df_matched, cases, controls,n_control_matched, n_case_matched, ratio, covs=[], group_col = 'Group', col_event = 'diag_vrs'):
    # Crear tabla de contingencia
    
    contingency_table = pd.crosstab(
        df_matched[col_event], #'diag_vrs' cama_bin
        df_matched['estado_inmunizacion'], 
        rownames=['Diagnóstico VRS'], 
        colnames=['Estado de Inmunización'], 
        margins=True
    )
    
    # Crear tabla de porcentajes
    percentage_table = contingency_table.div(contingency_table.loc['All'], axis=1) * 100
    percentage_by_column = contingency_table.div(contingency_table.sum(axis=0), axis=1) * 100
    
    if (0 not in contingency_table.columns) or (
        (contingency_table
        .drop(index='All', errors='ignore')
        .drop(columns='All', errors='ignore') == 0).any().any()):
        
        return 'No hay no inmunizados'
    
    # Extraer valores específicos para el cálculo del odds ratio
    a = contingency_table.loc[1, 1]  # Casos expuestos
    b = contingency_table.loc[1, 0]  # Casos no expuestos
    c = contingency_table.loc[0, 1]  # Controles expuestos
    d = contingency_table.loc[0, 0]  # Controles no expuestos

    
    # Calcular el odds ratio y la efectividad
    odds_ratio = (a * d) / (b * c)
    efectividad = (1 - odds_ratio) * 100

    def build_exog_condlogit(df_matched, covs, group_col):
        X_binary_main = df_matched[['estado_inmunizacion']].copy()

        if not covs:
            exog = X_binary_main
            return exog

        binarias = []
        numericas = []
        categoricas = []

        for col in covs:
            s = df_matched[col]

            # binaria (0/1) aunque venga con NaN
            vals = s.dropna().unique()
            if len(vals) > 0 and set(vals).issubset({0, 1}):
                binarias.append(col)
                continue

            # numérica
            if pd.api.types.is_numeric_dtype(s):
                numericas.append(col)
                continue

            # categórica (object/category/bool/otros)
            categoricas.append(col)

        X_bin = df_matched[binarias].copy() if binarias else pd.DataFrame(index=df_matched.index)

        # Escalar numéricas
        if numericas:
            scaler = StandardScaler()
            X_num_scaled = pd.DataFrame(
                scaler.fit_transform(df_matched[numericas]),
                columns=numericas,
                index=df_matched.index
            )
        else:
            X_num_scaled = pd.DataFrame(index=df_matched.index)

        # One-hot categóricas
        if categoricas:
            X_cat = pd.get_dummies(
                df_matched[categoricas].astype("category"),
                prefix=categoricas,
                drop_first=True,      # evita dummy trap
                dummy_na=False        # pon True si quieres una categoría explícita para NA
            )
        else:
            X_cat = pd.DataFrame(index=df_matched.index)

        exog = pd.concat([X_binary_main, X_bin, X_num_scaled, X_cat], axis=1)

        keep_cols = [
            c for c in exog.columns
            if df_matched.groupby(group_col)[c].nunique(dropna=False).gt(1).any()
        ]
        exog = exog[keep_cols]

        # (opcional) asegurar float
        exog = exog.astype(float)

        return exog

    # Uso en tu función:
    exog = build_exog_condlogit(df_matched, covs=covs, group_col=group_col)
    endog = df_matched[col_event]
    
    # X_binary_main = df_matched[['estado_inmunizacion']]

    # # Verificar si hay covariables
    # if covs:
    #     binarias = [col for col in covs if df_matched[col].dropna().isin([0, 1]).all()]
    #     no_binarias = [col for col in covs if col not in binarias]

    #     X_binarias = df_matched[binarias] if binarias else pd.DataFrame(index=df_matched.index)
    #     X_no_binarias = df_matched[no_binarias] if no_binarias else pd.DataFrame(index=df_matched.index)

    #     if not X_no_binarias.empty:
    #         scaler = StandardScaler()
    #         X_no_binarias_scaled = pd.DataFrame(scaler.fit_transform(X_no_binarias),
    #                                             columns=no_binarias, index=df_matched.index)
    #     else:
    #         X_no_binarias_scaled = pd.DataFrame(index=df_matched.index)
    #     exog = pd.concat([X_binary_main, X_binarias, X_no_binarias_scaled], axis=1)
        
    # else:
    #     exog = X_binary_main.copy()

    # endog = df_matched[col_event]

    model = ConditionalLogit(endog, exog, groups=df_matched[group_col])
    result = model.fit()
    
    # Crear DataFrame con coeficientes, odds ratios y efectividad
    odds_ratios = (
        pd.DataFrame({
            'Coeficientes': result.params[:1],
            'OR': np.exp(result.params[:1]),
            'Efectividad': (1 - np.exp(result.params[:1])) * 100,
        })
        .merge((1 - np.exp(result.conf_int()))[:1].rename({0:'0.975', 1: '0.025'}, axis=1)*100, left_index=True, right_index=True)
    )
    
    # Preparar los resultados en un formato adecuado para guardar en el diccionario
    results_df = pd.DataFrame({
        'Coeficientes': result.params[:1],
        'OR': np.exp(result.params[:1]),
        'Efectividad': (1 - np.exp(result.params[:1])) * 100,
        '0.975 Conf Interval': odds_ratios['0.975'],
        '0.025 Conf Interval': odds_ratios['0.025'],
        'Odds Ratio Manual': [odds_ratio],
        'Efectividad Manual': [efectividad]
    }).T  # 
    
    return {
        'df_matched': df_matched,
        'covariates': covs,
        'total_cases':cases.shape[0],
        'total_controls': controls.shape[0],
        'cases_matched': n_case_matched,
        'controls_matched': n_control_matched,
        'ratio': ratio,
        'contingency_table': contingency_table,
        'percentage_table': percentage_table,
        'percentage_by_column': percentage_by_column,
        'odds_ratio': odds_ratio,
        'efectividad': efectividad,
        'odds_ratios': odds_ratios,
        'results_df': results_df,
        'model_summary': result.summary(),
        'IC_distance': (result.conf_int()[1] - result.conf_int()[0]).values[0],
        'cases_no_inmune': b,
        'cases_inmune': a,
        'controles_no_inmune': d,
        'controles_inmune': c,
        'cases_column': f'{a}/{n_case_matched} ({(float(a)/float(n_case_matched))*100:.1f})',
        'controls_column': f'{c}/{n_control_matched} ({(float(c)/float(n_control_matched))*100:.1f})',
    }

def mylogit(df,cases, controls, ratio, n_control_matched, n_case_matched,prints=False, covs=[]):
    
    if prints:
        display(comparar_medias_test(df[~df.Matched_Case_RUN.isna()],df[df.Matched_Case_RUN.isna()],['PESO','SEMANAS','SEXO','edad_relativa']))
        print('\n')
        print("Dado tu estado de inmunizacion cual es probabilidad de estar contagiado")
        print(df.groupby('estado_inmunizacion')['diag_vrs'].mean())
        print('\n')
    #  analyze_vrs_data(df_matched, cases, controls,n_control_matched, n_case_matched, ratio, covs=[], group_col = 'Group', col_event = 'diag_vrs')   
    results = analyze_vrs_data(df_matched=df, 
                               cases=cases, 
                               controls=controls,
                               ratio=ratio,
                               n_control_matched=n_control_matched, 
                               n_case_matched=n_case_matched, 
                               covs=covs)
    if results == 'No hay no inmunizados':
        return results
    if prints:
        print("Tabla de Contingencia:\n", results['contingency_table'])
        print("\nTabla de Porcentajes:\n", results['percentage_table'])
        print(f"\nOdds Ratio: {results['odds_ratio']: .2}")
        print(f"Efectividad: {results['efectividad']: .3} %")
        print("\nResumen del modelo:\n", results['model_summary'])
    print("\nOdds Ratios y Efectividad:\n", results['odds_ratios'])
    print('\nIC disntace:', results['IC_distance'])
    
    return results

def charly_mip(df_cases, df_control, distance_vars, exact_var_match = ['sexo','region'], ratio="1:1"):
    
    r_1 = int(ratio.split(":")[0])
    r_2 = int(ratio.split(":")[1])
    
    cases = df_cases.copy().reset_index(drop=True)
    controls = df_control.copy().reset_index(drop=True)
    

    if 'fecha_nac' in distance_vars:
        ref_date = min(cases['fecha_nac'].min(), controls['fecha_nac'].min())

        cases['fecha_nac_dias'] = (cases['fecha_nac'] - ref_date).dt.days
        controls['fecha_nac_dias'] = (controls['fecha_nac'] - ref_date).dt.days

        # Actualizar distance_vars con la nueva columna numérica
        distance_vars = [var if var != 'fecha_nac' else 'fecha_nac_dias' for var in distance_vars]
    
    feasible_controls = {}
    star_time = time.time()
    for i, case in cases.iterrows():

        controls_copy = df_control.copy().reset_index(drop=True)
        
        for var in exact_var_match:
            controls_copy = controls_copy[controls_copy[var] == case[var]]
        
        feasible_controls[i] = controls_copy.index.tolist() if not controls_copy.empty else []
    end_time = time.time()
    
    print(f"creacion conjuntos A_i time: {end_time - star_time}")
    star_time = time.time()

    model = Model("Maximize_Matching")
    model.setParam('OutputFlag', 0)  
    x = {}
    d = {}
    y = {}
    inv_cov_matrix = np.linalg.inv(np.cov(controls[distance_vars].T))
    control_set = set()
    for i, control_indices in feasible_controls.items():
        
        y[i] = model.addVar(vtype=GRB.BINARY, name=f"y_{i}")
        
        for j in control_indices:
            d[(i, j)] = mahalanobis(cases.loc[i, distance_vars], controls.loc[j, distance_vars], inv_cov_matrix)
            x[(i, j)] = model.addVar(vtype=GRB.BINARY, name=f"x_{i}_{j}")   
            control_set.add(j)
    if r_1>1:
        z = {j: model.addVar(vtype=GRB.BINARY, name=f"z_{j}") for j in control_set}
         
    end_time = time.time()
    print(f"creacion variables time: {end_time - star_time}")
    distance_max = max(d.values())
    
    model.setObjective(sum((d[(i, j)] - distance_max)*x[(i, j)] for (i, j) in x.keys()), GRB.MINIMIZE)

    for i in feasible_controls:
        if r_2==1:
            model.addConstr(sum(x[(i, j)] for j in feasible_controls[i]) <= 1, name=f"max_controls_case_{i}")
        else:
            model.addConstr(sum(x[(i, j)] for j in feasible_controls[i]) == y[i]*r_2, name=f"max_controls_case_{i}")
        
        # for j in feasible_controls[i]:
        #     model.addConstr(x[(i, j)] <= y[i])
        #     model.addConstr(x[(i, j)] <= z[j])
           # model.addConstr(x[(i, j)] >= y[i] + z[j] - 1) #aporta???
            
    for j in controls_copy.index:
        if r_1==1:
            model.addConstr(sum(x[(i, j)] for i in feasible_controls if j in feasible_controls[i]) <= 1, name=f"max_case_control_{j}")
        else:
            model.addConstr(sum(x[(i, j)] for i in feasible_controls if j in feasible_controls[i]) == z[j]*r_1, name=f"max_case_control_{j}")

    star_time = time.time()
    model.optimize()
    
    end_time = time.time()
    print(f"optimize model time: {end_time - star_time}")
    
    cases_m =[]
    matched_pairs = []
    
    star_time = time.time()
    for (i, j), var in x.items():
        if var.X > 0.5:  # Solo emparejar si la variable de decisión es 1
                selected_control = controls.iloc[[j]].copy()
                selected_control['Matched_Case_RUN'] = cases.loc[i, 'RUN']
                matched_pairs.append(selected_control)
                cases_m.append(cases.loc[i, 'RUN'])
                
    end_time = time.time()
    print(f"matched_data time: {end_time - star_time}")    
            
    # Construir el DataFrame `matched_data` en el mismo formato
    matched_controls = pd.concat(matched_pairs, axis=0) if matched_pairs else pd.DataFrame()
    cases_matched = cases[cases['RUN'].isin(matched_controls['Matched_Case_RUN'])]
    matched_data = pd.concat([cases_matched, matched_controls], axis=0)
    matched_data['Group'] = matched_data['Matched_Case_RUN'].fillna(matched_data['RUN'])
    
    print(f'Total cases matched is : {len(cases_m)}')
    
    return matched_data, feasible_controls

def charly_double_mip(df_cases, df_control, distance_vars, exact_var_match = ['sexo','region'], ratio="1:1"):
    
    r_1 = int(ratio.split(":")[0])
    r_2 = int(ratio.split(":")[1])
    
    cases = df_cases.copy().reset_index(drop=True)
    controls = df_control.copy().reset_index(drop=True)
    

    if 'fecha_nac' in distance_vars:
        ref_date = min(cases['fecha_nac'].min(), controls['fecha_nac'].min())

        cases['fecha_nac_dias'] = (cases['fecha_nac'] - ref_date).dt.days
        controls['fecha_nac_dias'] = (controls['fecha_nac'] - ref_date).dt.days

        # Actualizar distance_vars con la nueva columna numérica
        distance_vars = [var if var != 'fecha_nac' else 'fecha_nac_dias' for var in distance_vars]
    
    feasible_controls = {}
    star_time = time.time()
    for i, case in cases.iterrows():

        controls_copy = df_control.copy().reset_index(drop=True)
        
        for var in exact_var_match:
            controls_copy = controls_copy[controls_copy[var] == case[var]]
        
        feasible_controls[i] = controls_copy.index.tolist() if not controls_copy.empty else []
    end_time = time.time()
    
    print(f"creacion conjuntos A_i time: {end_time - star_time}")
    
    model = Model("max_match")
    model.setParam('OutputFlag', 0)  
    x = {}
    d = {}
    y = {}
    
    if len(distance_vars) == 1:
        var = np.var(controls[distance_vars[0]], ddof=1)
        inv_cov_matrix = 1.0 / var if var > 0 else 0.0  # Protección por si varianza es 0
    else:
        inv_cov_matrix = np.linalg.inv(np.cov(controls[distance_vars].T))
        
    control_set = set()
    for i, control_indices in feasible_controls.items():
        y[i] = model.addVar(vtype=GRB.BINARY, name=f"y_{i}")
        
        for j in control_indices:
            if len(distance_vars) == 1:
                diff = cases.loc[i, distance_vars[0]] - controls.loc[j, distance_vars[0]]
                dist = np.sqrt(diff**2 * inv_cov_matrix)
            else:
                dist = mahalanobis(cases.loc[i, distance_vars], controls.loc[j, distance_vars], inv_cov_matrix)
            
            d[(i, j)] = dist
            x[(i, j)] = model.addVar(vtype=GRB.BINARY, name=f"x_{i}_{j}")
            control_set.add(j)

    if r_1 > 1:
        z = {j: model.addVar(vtype=GRB.BINARY, name=f"z_{j}") for j in control_set}

    end_time = time.time()
    print(f"creacion variables time: {end_time - star_time}")
    distance_max = max(d.values())
    
    print(distance_max)
    
    model.setObjective(sum(x[(i, j)] for (i, j) in x.keys()), GRB.MAXIMIZE)

    for i in feasible_controls.keys():
        if r_2==1:
            model.addConstr(sum(x[(i, j)] for j in feasible_controls[i]) <= 1, name=f"max_controls_case_{i}")
        else:
            model.addConstr(sum(x[(i, j)] for j in feasible_controls[i]) == y[i]*r_2, name=f"max_controls_case_{i}")
        
            
    for j in control_set:#control_set: #: #in controls_copy (asi estaba antes del 07-04)    ,   controls.index
        if r_1==1:
            model.addConstr(sum(x[(i, j)] for i in feasible_controls.keys() if j in feasible_controls[i]) <= 1, name=f"max_case_control_{j}")
        else:
            model.addConstr(sum(x[(i, j)] for i in feasible_controls.keys() if j in feasible_controls[i]) == z[j]*r_1, name=f"max_case_control_{j}")

    star_time = time.time()
    model.optimize()
    
    end_time = time.time()
    print(f"optimize model 1 time: {end_time - star_time}")
    
    optimal_value_model1 = model.objVal
    
    model.setObjective(sum((d[(i, j)] - distance_max)*x[(i, j)] for (i, j) in x.keys()), GRB.MINIMIZE)
    
    model.addConstr(sum(x[(i, j)] for (i, j) in x.keys()) >= optimal_value_model1, name="optim_m1")
    
    star_time = time.time()
    model.optimize()
    end_time = time.time()
    print(f"optimize model 2 time: {end_time - star_time}")
    
    star_time = time.time()
    
    cases_m =[]
    matched_pairs = []
    
    star_time = time.time()
    for (i, j), var in x.items():
        if r_1<=r_2:
            if var.X > 0.5:  # Solo emparejar si la variable de decisión es 1
                    selected_control = controls.iloc[[j]].copy()   ######## why not controls_copy??
                    selected_control['Matched_Case_RUN'] = cases.loc[i, 'RUN']
                    matched_pairs.append(selected_control)
                    cases_m.append(cases.loc[i, 'RUN'])
        # else:
        #     if var.X > 0.5:  # Solo emparejar si la variable de decisión es 1
        #             selected_control = controls.iloc[[i]].copy()   ######## why not controls_copy??
        #             selected_control['Matched_Case_RUN'] = cases.loc[j, 'RUN']
        #             matched_pairs.append(selected_control)
        #             cases_m.append(cases.loc[, 'RUN'])
                
    end_time = time.time()
    print(f"matched_data time: {end_time - star_time}")    
    
    
    matched_controls = pd.concat(matched_pairs, axis=0) if matched_pairs else pd.DataFrame()
    cases_matched = cases[cases['RUN'].isin(matched_controls['Matched_Case_RUN'])]
    matched_data = pd.concat([cases_matched, matched_controls], axis=0).reset_index(drop=True)
    matched_data['Group'] = matched_data['Matched_Case_RUN'].fillna(matched_data['RUN'])
    cases_m = cases_matched['RUN'].tolist()
    
    print(f'Total cases = {cases.shape[0]}, Total controls = {controls.shape[0]}')
    print(f'Total cases matched is : {len(cases_m)}, Total control matched is : {matched_controls.shape[0]}') 
    print('ratio: ' + ratio) 
    
    n_control_matched = matched_controls.shape[0]
    n_case_matched = len(cases_m)
    
    
    return matched_data, feasible_controls, cases, controls, ratio, n_control_matched, n_case_matched

def results_match(df_case_study,df_control_study,filtros_dic,match_vars_distance_nn,match_vars_exact_nn,match_vars_distance_IP,match_vars_exact_IP,weights, list_experiments = [],nn=False,mip=True, ratio='1:3',covs=[],dic_ratio={}):
    
    ruts_cariopaticos_1 = pd.read_csv(path_data/'lista_pacientes_riesgo.csv').query('card1')
    ruts_cariopaticos_2 = pd.read_csv(path_data/'lista_pacientes_riesgo.csv').query('card2')
    ruts_displasia = pd.read_csv(path_data/'lista_pacientes_riesgo.csv').query('displ')
    
    with open(path_data/'lista_ruts_cardio.pkl', 'rb') as f:
        lista_ruts_cardio = pickle.load(f)
    
    def variables_que_cambian(df, lista_vars):
        return [var for var in lista_vars if df[var].nunique(dropna=False) > 1]
       
    covs_copy = covs.copy()
    for name_filt, filtro in filtros_dic.items():
        # df = df_estudio.query(filtro).copy()
        # df_case = df.query("diag_vrs==True").copy()
        # df_control = df.query("diag_vrs==False").copy()
        if dic_ratio == {}:
            ratio = ratio
        else:
            ratio = dic_ratio[name_filt]
        
        r_1 = int(ratio.split(":")[0])
        r_2 = int(ratio.split(":")[1])
        
        if df_case_study.shape[0] < 3*df_control_study.shape[0]:
            control_group_name = 'all_born'
        else:
            control_group_name = 'francia'
        
        if r_1 > r_2:
            df_case = df_control_study.query("diag_vrs==False").query(filtro).copy()
            df_control = df_case_study.query("diag_vrs==True").query(filtro).copy()
        else:
            df_case = df_case_study.query("diag_vrs==True").query(filtro).copy()
            df_control = df_control_study.query("diag_vrs==False").query(filtro).copy()
        
        covs_check1 = variables_que_cambian(df_case, covs_copy)
        covs_check2 = variables_que_cambian(df_control, covs_copy)
        
        # if covs_check1 != covs_check2:
        #     print(covs_check1, covs_check2)
        #     raise ValueError("Las covariables no son las mismas en los dos dataframes")
        # else:
        
        if len(covs_check1) < len(covs_check2):
            covs = covs_check1.copy()
        else:
            covs = covs_check2.copy()
        
        n_casos = df_case.shape[0]
        n_controls = df_control.shape[0]

        print(name_filt, n_casos, n_controls, covs )

        if n_casos < 1 or n_controls < 1:
            continue
        if nn==True:
            
            
            lista_df= match_nn_max_dist_weigths(df_control, df_case,
                                                match_vars_nn= match_vars_distance_nn, 
                                                match_vars_exact = match_vars_exact_nn,
                                                match_vars_onehot=[],
                                                ratio=ratio,
                                                max_distance=5,
                                                weights=weights)
            if lista_df == "No hay controles disponibles para emparejar":
                dic_aux = {'filtro': name_filt,
                            'modelo': 'nn',
                            'problema': 'No hay controles'}
                list_experiments.append(dic_aux)
            else:
                matched_data, matched_incompleto, matched_controls, cases, controls, ratio, n_control_matched, n_case_matched = lista_df
            
                results = mylogit(matched_data,cases, controls, ratio, n_control_matched, n_case_matched, covs=covs)
                if results == 'No hay no inmunizados':
                    dic_aux = {'filtro': name_filt,
                            'modelo': 'nn',
                            'problema': 'No hay no inmunizados'}
                    list_experiments.append(dic_aux)
                else:
                    results['filtro'] = name_filt
                    results['model'] = 'nn'
                    list_experiments.append(results.copy())
        if mip==True:
            #try:
            matched_data, feasible_controls, cases, controls, ratio, n_control_matched, n_case_matched = charly_double_mip(df_cases=df_case,
                                                            df_control=df_control, 
                                                            distance_vars=match_vars_distance_IP, 
                                                            exact_var_match = match_vars_exact_IP, 
                                                            ratio=ratio)
            name_base = 'matched_data' + '_' + name_filt + '_' + control_group_name + '.csv'
            matched_data.to_csv(path_data/name_base)
            results = mylogit(df=matched_data,
                              cases=cases, 
                              controls=controls, 
                              ratio=ratio, 
                              n_control_matched=n_control_matched, 
                              n_case_matched=n_case_matched, 
                              covs=covs)
            if results == 'No hay no inmunizados':
                dic_aux = {'filtro': name_filt,
                        'modelo': 'MIP',
                        'problema': 'No hay no inmunizados'}
                list_experiments.append(dic_aux)
            else:
                results['filtro'] = name_filt
                results['model'] = 'MIP'
                list_experiments.append(results.copy())
            # except:
            #     print('No se pudo hacer el matching con mip se hizo nn')
            #     lista_df= match_nn_max_dist_weigths(df_control, df_case,
            #                                         match_vars_nn= match_vars_distance_nn, 
            #                                         match_vars_exact = match_vars_exact_nn,
            #                                         match_vars_onehot=[],
            #                                         ratio=ratio,
            #                                         max_distance=100,
            #                                         weights=weights)
            #     if lista_df == "No hay controles disponibles para emparejar":
            #         dic_aux = {'filtro': name_filt,
            #                 'modelo': 'nn',
            #                 'problema': 'No hay controles'}
            #         list_experiments.append(dic_aux)
            #     else:
            #         matched_data, matched_incompleto, matched_controls, cases, controls, ratio, n_control_matched, n_case_matched = lista_df

            #         results = mylogit(matched_data,cases, controls, ratio, n_control_matched, n_case_matched, covs=covs)
            #         if results == 'No hay no inmunizados':
            #             dic_aux = {'filtro': name_filt,
            #                     'modelo': 'nn',
            #                     'problema': 'No hay no inmunizados'}
            #             list_experiments.append(dic_aux)
            #         else:
            #             results['filtro'] = name_filt
            #             results['model'] = 'nn'
            #             list_experiments.append(results.copy())
    return list_experiments

def tabla_marcel(list_experiments):
    
    experiment_dfs = []
    list_experiments_edit = [dic for dic in list_experiments if len(dic) >4] #!= 3

    for result in list_experiments_edit:
        results_df_transposed =(
            result['results_df'].T
            .assign(
                efectividad = lambda df : np.where(df["0.025 Conf Interval"].values[0].round(2) < -100, 
                                                'na',
                                                f'{str(df["Efectividad"].values[0].round(2))} ({str(df["0.025 Conf Interval"].values[0].round(2))}; {str(df["0.975 Conf Interval"].values[0].round(2))})'),    
                OR_w_IC = lambda df : f'{str(df["OR"].values[0].round(2))} ({str((100 - df["0.975 Conf Interval"].values[0].round(2)).round(2))}; {str((100 - df["0.025 Conf Interval"].values[0].round(2)).round(2))})',
                        
            )
            .assign(
                    total_cases = result['total_cases'],
                    inmune_cases = result['cases_inmune'],
                    total_controls = result['total_controls'],
                    inmune_controls = result['controles_inmune'],
                    cases_matched = result['cases_matched'],
                    controls_matched = result['controls_matched'],
                    cases_column = result['cases_column'],
                    controls_column = result['controls_column'],
                    ratio = result['ratio'],
                    method = result['model'],
                    filtro = result['filtro']
                )
            .pipe(lambda df:df.set_index(['filtro'])) 
            ) 
        
        experiment_dfs.append(results_df_transposed)

    all_results_df = pd.concat(experiment_dfs, ignore_index=False)

    marcel_summary= (all_results_df[['efectividad','OR_w_IC']].reset_index())
    
    return marcel_summary, all_results_df

def tabla_final(all_results_df,marcel_summary):
    compare = all_results_df[['total_cases','inmune_cases', 'total_controls','inmune_controls','cases_matched', 'controls_matched','ratio','cases_column','controls_column']]
    # compare.columns = ['total_cases','inmune_cases', 'total_controls','inmune_controls','controls_matched', 'ratio','cases_matched','cases_column','controls_column']
    # compare = compare.assign(prop_cases_match = lambda x: x.cases_matched/x.total_cases,
                            # prop_controls_match = lambda x: x.controls_matched/x.total_controls)

    # compare = compare[['total_cases','cases_matched','inmune_cases','controls_matched','inmune_controls','cases_column','controls_column']] #['total_cases','cases_matched', 'total_controls','controls_matched']
    compare = compare[['total_cases','cases_column','controls_column']] #['total_cases','cases_matched', 'total_controls','controls_matched']
    enes_matches = compare.rename(columns={'cases_column':'Cases_immune/matched','controls_column':'Controls_immune/matched'}).round(2).reset_index() #.unstack(-1).round(2)
    # all_n_eff = marcel_summary[['filtro','efectividad']].merge(enes_matches, how='left', on='filtro')
    all_n_eff = enes_matches.merge(marcel_summary[['filtro','OR_w_IC','efectividad']], how='left', on='filtro')
    return enes_matches, all_n_eff

def summary_eff(list_experiments_all_born, list_experiments_francia_wPrevi, list_experiments_francia_not_previ):

    marcel_summary_all_born, all_results_df_all_born = tabla_marcel(list_experiments_all_born)
    enes_matches_all_born, all_n_eff_all_born  = tabla_final(all_results_df_all_born,marcel_summary_all_born)
    df_all_born = all_n_eff_all_born

    marcel_summary_francia_wPrevi, all_results_df_francia_wPrevi = tabla_marcel(list_experiments_francia_wPrevi)
    enes_matches_francia_wPrevi, all_n_eff_francia_wPrevi  = tabla_final(all_results_df_francia_wPrevi,marcel_summary_francia_wPrevi)
    df_francia_previ = all_n_eff_francia_wPrevi

    marcel_summary_francia_not_previ, all_results_df_francia_not_previ = tabla_marcel(list_experiments_francia_not_previ)
    enes_matches_francia_not_previ, all_n_eff_francia_not_previ  = tabla_final(all_results_df_francia_not_previ,marcel_summary_francia_not_previ)
    df_francia_not_previ = all_n_eff_francia_not_previ

    # Paso 1: Asegúrate de que el índice sea 'filtro' para todas
    df_all_born_macro = df_all_born.set_index("filtro")
    df_francia_previ_macro = df_francia_previ.set_index("filtro")
    df_francia_not_previ_macro = df_francia_not_previ.set_index("filtro")

    # Paso 2: Renombrar columnas con MultiIndex
    df_all_born_macro.columns = pd.MultiIndex.from_product([["all_born"], df_all_born_macro.columns])
    df_francia_previ_macro.columns = pd.MultiIndex.from_product([["francia_previ"], df_francia_previ_macro.columns])
    df_francia_not_previ_macro.columns = pd.MultiIndex.from_product([["francia_not_previ"], df_francia_not_previ_macro.columns])

    # Paso 3: Concatenar horizontalmente por el índice (filtro)
    df_final = pd.concat([df_all_born_macro, df_francia_not_previ_macro, df_francia_previ_macro], axis=1)
    # df_final.columns = pd.MultiIndex.from_tuples([
    #     (nivel_superior, "prop" if nivel_inferior == "prop_cases_match" else nivel_inferior)
    #     for nivel_superior, nivel_inferior in df_final.columns
    # ])
    df_final.columns = pd.MultiIndex.from_tuples([
        (nivel_superior, "matched" if nivel_inferior == "cases_matched" else nivel_inferior)
        for nivel_superior, nivel_inferior in df_final.columns
    ])
    df_final.columns = pd.MultiIndex.from_tuples([
        (nivel_superior, "cases" if nivel_inferior == "total_cases" else nivel_inferior)
        for nivel_superior, nivel_inferior in df_final.columns
    ])
    return df_final

def tabla_1_match(df, id_col='RUN', group_col='Group', case_label='Caso', control_label='Control', decimales=2):

    df = df.copy()
    df['grupo_tabla'] = np.where(df[id_col] == df[group_col], case_label, control_label)
    df['sexo'] = df['sexo'].map({0: 'Female', 1: 'Male'})

    variables_numericas_mediana = ['edad_relativa']
    variables_numericas_media = ['SEMANAS', 'PESO']
    variables_categoricas = ['sexo']
    macrozona_var = 'Macrozones'

    etiquetas_variables = {
        'edad_relativa': 'Age (days)',
        'SEMANAS': 'Gestational age at birth (weeks)',
        'PESO': 'Birth weight (grams)',
        'sexo': 'Sex',
        'Macrozones': 'Macro-zone',
        'cardio1': 'Congenital heart disease',
        'prematuro': 'Prematurity',
        'event_upc': 'PICU admission'
    }

    grupos = df['grupo_tabla'].unique()
    tabla = []
    index = []
    n_por_grupo = {g: len(df[df['grupo_tabla'] == g]) for g in grupos}
    nombres_columnas = {g: f"{g} (n = {n_por_grupo[g]})" for g in grupos}

    for var in variables_numericas_mediana:
        fila = {}
        for g in grupos:
            subset = df[df['grupo_tabla'] == g][var]
            mediana = subset.median()
            iqr = subset.quantile(0.75) - subset.quantile(0.25)
            fila[nombres_columnas[g]] = f'{mediana:.{decimales}f} ({iqr:.{decimales}f})'
        index.append((etiquetas_variables.get(var, var), ''))
        tabla.append(fila)

    for var in variables_categoricas:
        index.append((etiquetas_variables.get(var, var), ''))
        tabla.append({})  # Fila vacía
        niveles = df[var].dropna().unique()
        for val in niveles:
            fila = {}
            for g in grupos:
                subset = df[df['grupo_tabla'] == g]
                total = len(subset)
                cuenta = (subset[var] == val).sum()
                porcentaje = 100 * cuenta / total if total > 0 else 0
                fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
            index.append((etiquetas_variables.get(var, var), val))
            tabla.append(fila)

    for var in variables_numericas_media:
        fila = {}
        for g in grupos:
            subset = df[df['grupo_tabla'] == g][var]
            media = subset.mean()
            sd = subset.std()
            fila[nombres_columnas[g]] = f'{media:.{decimales}f} ({sd:.{decimales}f})'
        index.append((etiquetas_variables.get(var, var), ''))
        tabla.append(fila)

    index.append(('Risk groups', ''))
    tabla.append({})
    for var in ['cardio1', 'prematuro']:
        fila = {}
        for g in grupos:
            subset = df[df['grupo_tabla'] == g]
            total = len(subset)
            cuenta = (subset[var] == 1).sum()
            porcentaje = 100 * cuenta / total if total > 0 else 0
            fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
        index.append(('Risk groups', etiquetas_variables.get(var, var)))
        tabla.append(fila)

    index.append((etiquetas_variables.get(macrozona_var, macrozona_var), ''))
    tabla.append({})
    niveles = df[macrozona_var].dropna().unique()
    for val in niveles:
        fila = {}
        for g in grupos:
            subset = df[df['grupo_tabla'] == g]
            total = len(subset)
            cuenta = (subset[macrozona_var] == val).sum()
            porcentaje = 100 * cuenta / total if total > 0 else 0
            fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
        index.append((etiquetas_variables.get(macrozona_var, macrozona_var), val))
        tabla.append(fila)

    fila = {}
    for g in grupos:
        subset = df[df['grupo_tabla'] == g]
        total = len(subset)
        cuenta = (subset['cama'] == 'UPC').sum()
        porcentaje = 100 * cuenta / total if total > 0 else 0
        fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
    index.append((etiquetas_variables.get('cama', 'PICU'), ''))
    tabla.append(fila)

    diag_labels = {
        'A090': 'Infectious diarrhoea',
        'R681': 'Other general symptoms and signs',
        'R11X': 'Nausea and vomiting',
        'R104': 'Abdominal pain',
        'S099': 'Other head injuries',
        'R634': 'Feeding difficulties',
        'R633': 'Polydipsia',
        'N390': 'Urinary tract infection',
        'P599': 'Other neonatal jaundice',
    }
    
    diag_agrupao = {
        'A099': 'Other gastroenteritis and colitis',
        'A099': 'Other gastroenteritis and colitis',
        'A080': 'Other gastroenteritis and colitis',
        'A083': 'Other gastroenteritis and colitis',
        'A082': 'Other gastroenteritis and colitis',
        'A084': 'Other gastroenteritis and colitis',
        'A090': 'Other gastroenteritis and colitis',
        'A081': 'Other gastroenteritis and colitis',
        'A085': 'Other gastroenteritis and colitis',
    }

    diag_frecuencias = []
    for diag_code, label in diag_labels.items():
        fila = {}
        count_control = 0
        for g in grupos:
            if g == case_label:
                fila[nombres_columnas[g]] = '-'  # No diagnóstico para los casos
            else:
                subset = df[(df['grupo_tabla'] == g) & (df['DIAG1'] == diag_code)]
                total = len(df[df['grupo_tabla'] == g])
                cuenta = len(subset)
                count_control = cuenta
                porcentaje = 100 * cuenta / total if total > 0 else 0
                fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
        diag_frecuencias.append((count_control, ('Reason for hospital visit (controls)', label), fila))
    
    
    count_control_agrupao = 0
    for diag_code, label in diag_agrupao.items():
        
        for g in grupos:
            if g == control_label:
                subset = df[(df['grupo_tabla'] == g) & (df['DIAG1'] == diag_code)]
                total = len(df[df['grupo_tabla'] == g])
                cuenta = len(subset)
                count_control_agrupao += cuenta
    fila_agrupao = {}  
    porcentaje_agrupao = 100 * count_control_agrupao / total if total > 0 else 0
    fila_agrupao[nombres_columnas[control_label]] = f'{count_control_agrupao} ({porcentaje_agrupao:.1f}%)'
    fila_agrupao[nombres_columnas[case_label]] = '-'
    diag_frecuencias.append((count_control_agrupao, ('Reason for hospital visit (controls)', list(diag_agrupao.values())[0]), fila_agrupao))
    
    index.append(('Reason for hospital visit (controls)', ''))
    tabla.append({})
    for _, idx, fila in sorted(diag_frecuencias, key=lambda x: x[0], reverse=True):
        index.append(idx)
        tabla.append(fila)

    multiindex = pd.MultiIndex.from_tuples(index, names=['Variable', 'Value'])
    return pd.DataFrame(tabla, index=multiindex)

def summary_eff_aux(
    list_experiments_all_born=None,
    list_experiments_francia_wPrevi=None,
    list_experiments_francia_not_previ=None
):
    import pandas as pd

    # Función auxiliar para procesar un set de experimentos y devolver un DataFrame formateado
    def process_experiments(list_experiments, label):
        # 1. Llamamos a las funciones tabla_marcel y tabla_final
        marcel_summary, all_results_df = tabla_marcel(list_experiments)
        enes_matches, all_n_eff = tabla_final(all_results_df, marcel_summary)

        # 2. Ajustamos índices y renombramos columnas
        df_macro = all_n_eff.set_index("filtro")
        df_macro.columns = pd.MultiIndex.from_product([[label], df_macro.columns])
        return df_macro

    # Acumulamos aquí los DataFrames procesados
    frames = []

    # Procesamos all_born si no es None
    if list_experiments_all_born is not None:
        df_all_born_macro = process_experiments(list_experiments_all_born, "all_born")
        frames.append(df_all_born_macro)

    # Procesamos francia_not_previ si no es None
    if list_experiments_francia_not_previ is not None:
        df_francia_not_previ_macro = process_experiments(list_experiments_francia_not_previ, "francia_not_previ")
        frames.append(df_francia_not_previ_macro)

    # Procesamos francia_wPrevi si no es None
    if list_experiments_francia_wPrevi is not None:
        df_francia_previ_macro = process_experiments(list_experiments_francia_wPrevi, "francia_previ")
        frames.append(df_francia_previ_macro)

    # Si no pasaste nada, devolvemos None o un aviso
    if not frames:
        print("No se recibió ninguna lista de experimentos")
        return None

    # Concatenamos horizontalmente
    df_final = pd.concat(frames, axis=1)

    # Renombramos las columnas finales según tu lógica
    # 1) 'cases_matched' -> 'matched'
    df_final.columns = pd.MultiIndex.from_tuples([
        (lvl1, "matched" if lvl2 == "cases_matched" else lvl2)
        for (lvl1, lvl2) in df_final.columns
    ])
    # 2) 'total_cases' -> 'cases'
    df_final.columns = pd.MultiIndex.from_tuples([
        (lvl1, "cases" if lvl2 == "total_cases" else lvl2)
        for (lvl1, lvl2) in df_final.columns
    ])

    return df_final

def df_problemas_matching(list_1=[],list_2=[],list_3=[]):

    copia_born = list_1.copy()
    copia_france_previ= list_2.copy()
    copia_france_not_previ = list_3.copy()

    # Filtrar los experimentos que tienen exactamente 3 elementos
    lista_problemas_AB = [dic for dic in copia_born if len(dic) == 3]
    lista_problemas_francia_previ = [dic for dic in copia_france_previ if len(dic) == 3]
    lista_problemas_francia_not_previ = [dic for dic in copia_france_not_previ if len(dic) == 3]

    # Agregar la columna 'controles' a cada lista de problemas
    for d in lista_problemas_AB:
        d['controles'] = 'all_born'

    for d in lista_problemas_francia_previ:
        d['controles'] = 'francia_previ'

    for d in lista_problemas_francia_not_previ:
        d['controles'] = 'francia_not_previ'
        
    # Crear el DataFrame final con todos los problemas
    df_problemas = pd.DataFrame(lista_problemas_AB + lista_problemas_francia_previ + lista_problemas_francia_not_previ)
    
    return df_problemas

def match_ps_max_dist_weights(
    df_control,
    df_case,
    match_vars_nn,
    match_vars_exact=[],
    match_vars_onehot=[],
    ratio="1:1",
    max_distance=0.02,          # caliper en PS (recomendado 0.01–0.05)
    weights=None,               # opcional: ponderar variables en el PS (multiplicando columnas)
    treat_label_case=0,         # por defecto df_control=treated(1) y df_case=0
):
    """
    Matching basado en propensity score (PS) con restricciones exactas y sin reemplazo.
    Devuelve los mismos argumentos que match_nn_max_dist_weigths.

    Nota: max_distance ahora es caliper de |ps_case - ps_control|.
    """
    if weights is None:
        weights = {}

    controls = df_control.copy().reset_index(drop=True)
    cases = df_case.copy().reset_index(drop=True)

    case_count, control_count = map(int, ratio.split(":"))
    matched_pairs = []
    matched_incompleto = []
    used_indices = set()

    # -------------------------------
    # Prints iniciales
    # -------------------------------
    print("========================================")
    print(f"[PS MATCH] ratio={ratio}  (case:control = {case_count}:{control_count})")
    print(f"[PS MATCH] Total cases    = {cases.shape[0]}")
    print(f"[PS MATCH] Total controls = {controls.shape[0]}")
    print("========================================")

    # ---------- 1) Construir matriz X para PS ----------
    # numéricas
    Xc_num = controls[match_vars_nn].copy() if match_vars_nn else pd.DataFrame(index=controls.index)
    Xk_num = cases[match_vars_nn].copy()    if match_vars_nn else pd.DataFrame(index=cases.index)

    # aplicar pesos a numéricas
    for col, w in weights.items():
        if col in Xc_num.columns:
            Xc_num[col] = Xc_num[col] * w
            Xk_num[col] = Xk_num[col] * w

    # categóricas a one-hot (fit conjunto)
    if match_vars_onehot:
        enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
        combined_cat = pd.concat([cases[match_vars_onehot], controls[match_vars_onehot]], axis=0)

        X_cat = enc.fit_transform(combined_cat.astype(str))
        cat_cols = enc.get_feature_names_out(match_vars_onehot)

        X_cat = pd.DataFrame(X_cat, columns=cat_cols, index=combined_cat.index)

        Xk_cat = X_cat.iloc[:len(cases)].reset_index(drop=True)
        Xc_cat = X_cat.iloc[len(cases):].reset_index(drop=True)
    else:
        Xk_cat = pd.DataFrame(index=cases.index)
        Xc_cat = pd.DataFrame(index=controls.index)

    # juntar X
    X_case = pd.concat([Xk_num.reset_index(drop=True), Xk_cat.reset_index(drop=True)], axis=1)
    X_ctrl = pd.concat([Xc_num.reset_index(drop=True), Xc_cat.reset_index(drop=True)], axis=1)

    X_all = pd.concat([X_case, X_ctrl], axis=0).fillna(0)

    # y: treated=1 para controls (por defecto), treated=0 para cases
    y_all = np.r_[np.zeros(len(cases), dtype=int), np.ones(len(controls), dtype=int)]
    if treat_label_case == 1:
        y_all = 1 - y_all

    # ---------- 2) Estimar PS ----------
    ps_model = LogisticRegression(max_iter=2000, n_jobs=-1, solver="lbfgs")
    ps_model.fit(X_all, y_all)
    ps_all = ps_model.predict_proba(X_all)[:, 1]

    cases["ps"] = ps_all[:len(cases)]
    controls["ps"] = ps_all[len(cases):]

    # -------------------------------
    # Prints de PS (distribución)
    # -------------------------------
    print("[PS MATCH] PS summary - cases:")
    print(cases["ps"].describe(percentiles=[0.05, 0.5, 0.95]).to_string())
    print("[PS MATCH] PS summary - controls:")
    print(controls["ps"].describe(percentiles=[0.05, 0.5, 0.95]).to_string())
    print("========================================")

    # ---------- helpers ----------
    def filter_exact(df, row, match_vars_exact):
        for var in match_vars_exact:
            df = df[df[var] == row[var]]
        return df

    # ---------- 3) Matching ----------
    if case_count <= control_count:
        neighbors_count = control_count

        for i, case_row in cases.iterrows():
            potential_controls = filter_exact(controls, case_row, match_vars_exact)
            if potential_controls.empty:
                matched_incompleto.append(case_row["RUN"])
                continue

            ps_pool = potential_controls[["ps"]].values
            nn = NearestNeighbors(n_neighbors=min(neighbors_count * 20, len(ps_pool)))
            nn.fit(ps_pool)

            distances, indices = nn.kneighbors([[case_row["ps"]]])

            selected_controls = []
            for idx_local, dist in zip(indices[0], distances[0]):
                ctrl_global_idx = potential_controls.index[idx_local]
                if dist <= max_distance and ctrl_global_idx not in used_indices:
                    selected_controls.append(controls.loc[[ctrl_global_idx]].copy())
                    used_indices.add(ctrl_global_idx)
                if len(selected_controls) == neighbors_count:
                    break

            if len(selected_controls) == neighbors_count:
                for control in selected_controls:
                    control["Matched_Case_RUN"] = case_row["RUN"]
                    matched_pairs.append(control)
            else:
                matched_incompleto.append(case_row["RUN"])

        if matched_pairs == []:
            print("[PS MATCH] No matched pairs found.")
            return "No hay controles disponibles para emparejar"

        matched_controls = pd.concat(matched_pairs, axis=0) if matched_pairs else pd.DataFrame()
        cases_matched = cases[cases["RUN"].isin(matched_controls["Matched_Case_RUN"])]

        matched_data = pd.concat([cases_matched, matched_controls], axis=0).reset_index(drop=True)
        matched_data["Group"] = matched_data["Matched_Case_RUN"].fillna(matched_data["RUN"])

        n_case_matched = cases_matched["RUN"].nunique()
        n_control_matched = matched_controls.shape[0]

    else:
        neighbors_count = case_count

        for j, ctrl_row in controls.iterrows():
            potential_cases = filter_exact(cases, ctrl_row, match_vars_exact)
            if potential_cases.empty:
                matched_incompleto.append(ctrl_row["RUN"])
                continue

            ps_pool = potential_cases[["ps"]].values
            nn = NearestNeighbors(n_neighbors=min(neighbors_count * 20, len(ps_pool)))
            nn.fit(ps_pool)

            distances, indices = nn.kneighbors([[ctrl_row["ps"]]])

            selected_cases = []
            for idx_local, dist in zip(indices[0], distances[0]):
                case_global_idx = potential_cases.index[idx_local]
                if dist <= max_distance and case_global_idx not in used_indices:
                    selected_cases.append(cases.loc[[case_global_idx]].copy())
                    used_indices.add(case_global_idx)
                if len(selected_cases) == neighbors_count:
                    break

            if len(selected_cases) == neighbors_count:
                for case in selected_cases:
                    case["Matched_Case_RUN"] = ctrl_row["RUN"]
                    matched_pairs.append(case)
            else:
                matched_incompleto.append(ctrl_row["RUN"])

        matched_controls = pd.concat(matched_pairs, axis=0) if matched_pairs else pd.DataFrame()
        cases_matched = controls[controls["RUN"].isin(matched_controls["Matched_Case_RUN"])]

        matched_data = pd.concat([cases_matched, matched_controls], axis=0).reset_index(drop=True)
        matched_data["Group"] = matched_data["Matched_Case_RUN"].fillna(matched_data["RUN"])

        n_control_matched = cases_matched["RUN"].nunique()
        n_case_matched = matched_controls.shape[0]

    # -------------------------------
    # Prints finales
    # -------------------------------
    print("========================================")
    print(f"[PS MATCH] Matched cases (unique RUN)   = {n_case_matched}")
    print(f"[PS MATCH] Matched controls (rows)      = {n_control_matched}")
    print(f"[PS MATCH] Unmatched primary units      = {len(matched_incompleto)}")
    print("========================================")

    return matched_data, matched_incompleto, matched_controls, cases, controls, ratio, n_control_matched, n_case_matched

