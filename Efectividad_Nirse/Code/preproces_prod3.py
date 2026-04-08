import pandas as pd
import pickle
from datetime import datetime, timedelta
from difflib import SequenceMatcher
import matplotlib.pyplot as plt
import os
from pathlib import Path
import seaborn as sns
import numpy as np
from lifelines.utils import to_long_format
from lifelines.utils import add_covariate_to_timeline
from lifelines import CoxTimeVaryingFitter
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import ks_2samp
from scipy import stats
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as mtick
import warnings
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from IPython.display import display
warnings.filterwarnings("ignore")

path_actual = Path.cwd()
path_data = path_actual.parent / 'Data'

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

dias = list(range(7, (40-13)*7, 7)) 
meses = [1,2,3,4,5,6]
diagnosticosVRS = ['J121', 'J205','J219', 'J210', 'B974']##['J120','J122','J123','J204','J206','J207','J211'] #['J121', 'J205','J219', 'J210', 'B974'] #['A099', 'A080', 'A083', 'A082', 'A084', 'A090', 'A081', 'A085'] pachecos_code, leos_code ['J120','J122','J123','J204','J206','J207','J211'] 
diagnosticos_upc = [406, 412, 415, 405, 411, 414, 310, 311, 312, 320, 323, 324]
months_list = ['Abril', 'Mayo', 'Junio', 'Julio', 'Agosto', 'Septiembre', 'Octubre']

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
              
              fecha_nac = lambda x: pd.to_datetime(x['FECHA_NACIMIENTO'], format='%d%b%Y'),
              fechaIng_any = lambda x: pd.to_datetime(x['FECHA_INGRESO'], format='%d%b%Y'),
              fechaEgr = lambda x: pd.to_datetime(x['FECHA_EGRESO'], format='%d%b%Y'),
              FECHA_INMUNIZACION = lambda x: pd.to_datetime(x['FECHA_INMUNIZACION'], format='%d%b%Y'),
              fechaInm = lambda x: pd.to_datetime(x['FECHA_INMUNIZACION'], format='%d%b%Y'),
              
              
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
              
              
              #features birth baby
              group = lambda x: np.where(x['fecha_nac'] < pd.to_datetime("2024-04-01"), "CATCH_UP", "SEASONAL"),
              sexo = lambda x: np.where(x['SEXO']==1.0,1,0),
              prematuro_extremo = lambda x: np.where((x['SEMANAS']>=22) & (x['SEMANAS']<=27),1,0),
              muy_prematuro = lambda x: np.where((x['SEMANAS']>=28) & (x['SEMANAS']<=32),1,0),
              prematuro_moderado = lambda x: np.where((x['SEMANAS']>=33) & (x['SEMANAS']<=36),1,0),
              prematuro = lambda x: np.where((x['SEMANAS']<=35),1,0),  ### CAMBIE ESTA WEA 15-05-2025 de 36 a 35
              atypic_mom_age = lambda x: np.where((x.EDAD_M>=45) | (x.EDAD_M<=19), 1, 0),
              region = lambda x: x['NOMBRE_REGION'].fillna('DESCONOCIDO').apply(lambda x: mapear_region(x, regiones) if isinstance(x, str) else x),
              Macrozona2 = lambda x: x['region'].map(region_a_macrozona2),
              
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
              vrs_pre_campaña = lambda x: np.where(x.fechaIng_vrs < pd.to_datetime("2024-04-01"), 1, 0), ################################### '2023-11-01'
              lrti_pre_campaña = lambda x: np.where(x.fechaIng_LRTI < pd.to_datetime("2024-04-01"), 1, 0),
              any_pre_campaña = lambda x: np.where(x.fechaIng_any < pd.to_datetime("2024-04-01"), 1, 0),
              upc_pre_campaña = lambda x: np.where(x.fecha_upc_vrs < pd.to_datetime("2024-04-01"), 1, 0)
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

def filtros_IH_new(df,dias=dias, meses=meses, semana_lower=24,semana_upper=42):

    # rut_eliminar =  ["bed99009d64eb031ead9235037fc95761d6f334e1d0bc27be4349f1734ca5b2f"]
    # df = df[~df.RUN.isin(rut_eliminar)]
    
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
    
    ################################################### ANALISIS MARCA ###############################
    
    # print('is_hapenin')
    # media_inmune_catch = np.mean((df1.query('group=="CATCH_UP"').fechaInm - df1.query('group=="CATCH_UP"').fecha_nac).dt.days)
    # media_inmune_nb = np.mean((df1.query('group=="SEASONAL"').fechaInm - df1.query('group=="SEASONAL"').fecha_nac).dt.days)
    
    # df1['random_mark'] = df1['MARCA']

    # marca_1_indices_catchup = df1.query('(group=="CATCH_UP") & (MARCA == 1)').RUN

    # num_to_change = int(0.111 * len(marca_1_indices_catchup))

    # random_run = np.random.choice(marca_1_indices_catchup, size=num_to_change, replace=False)


    # df1.loc[df1.RUN.isin(random_run), 'random_mark'] = 0
    
    # marca_1_indices_nb = df1.query('(group=="SEASONAL") & (MARCA == 1)').RUN

    # num_to_change_nb = int(0.053 * len(marca_1_indices_nb))

    # random_run_nb = np.random.choice(marca_1_indices_nb, size=num_to_change_nb, replace=False)

    # # Cambiar estos índices a 0 en la columna 'random_mark'
    # df1.loc[df1.RUN.isin(random_run_nb), 'random_mark'] = 0
    
    # df1['fechaInm'] = np.where((df1.random_mark==1) & (df1.group=="SEASONAL"), df1.fecha_nac + pd.DateOffset(days = media_inmune_nb), np.where((df1.random_mark==1) & (df1.group=="CATCH_UP"), df1.fecha_nac + pd.DateOffset(days = media_inmune_catch), df1.fechaInm))
    
    ##################################################################################
    
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
            df_onedit = df_onedit[~df_onedit.RUN.isin(df_onedit[(df_onedit['fechaIng_vrs'] < pd.to_datetime("2025-03-01")) & 
                                                                (df_onedit['fechaIng_vrs'] > pd.to_datetime("2024-09-30")) &
                                                                (df_onedit['fecha_nac'] > pd.to_datetime("2024-09-30"))].RUN.unique())]
        else:
            df_onedit = df_onedit[~df_onedit.RUN.isin(df_onedit[df_onedit[key] < pd.to_datetime("2024-04-01")].RUN.unique())]
            df_onedit = df_onedit[~df_onedit.RUN.isin(df_onedit[(df_onedit[key] < pd.to_datetime("2025-03-01")) & 
                                                                (df_onedit[key] > pd.to_datetime("2024-09-30")) &
                                                                (df_onedit['fecha_nac'] > pd.to_datetime("2024-09-30"))].RUN.unique())]
        
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
    
    fecha_dt = pd.to_datetime('2025-09-30') #+ pd.DateOffset(days=7) 2024-09-30
    
    for key, df_f in dfs_fecha.items():
    
        for m in meses:
            nombre_columna = f'age_{m}m'
            nombre_columna_bin = f'si_{m}_meses'
            df_f[nombre_columna] = df_f['fecha_nac'] + pd.DateOffset(months=m)
            df_f.loc[(df_f[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
            df_f[nombre_columna_bin] = 1 - df_f[nombre_columna].isna()
        
        
        for d in dias:
            nombre_columna = f'inm_{d}d'
            nombre_columna_bin = f'inm_mayor_{d}d'
            df_f[nombre_columna] = df_f['fechaInm'] + pd.DateOffset(days=d)
            df_f.loc[(df_f[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
            df_f[nombre_columna_bin] = 1 - df_f[nombre_columna].isna()
        

        for fecha_col, event_col in fechas_events.items():
            df_f.loc[df_f[fecha_col] > fecha_dt, fecha_col] = pd.NaT
            df_f[event_col] = df_f[fecha_col].notnull().astype(int)
    
        ##########################################3 ESTO LO CAMBIÉ EL 26 DE MAYO DE 2025 ##########################################################
        
       # df_save = df_f.loc[(df_f.fecha_nac <= pd.to_datetime("2024-09-30")) | (df_f.fecha_nac.dt.isocalendar().year == 2023)] 
        df_save = df_f.copy()
        ############################################################################################################################################
        dfs_fecha_actualiza[key] = df_save
        
    df_filtrado_any = dfs_fecha_actualiza['fechaIng_any']
    df_filtrado_LRTI = dfs_fecha_actualiza['fechaIng_LRTI']
    df_filtrado_vrs = dfs_fecha_actualiza['fechaIng_vrs']
    df_filtrado_upc = dfs_fecha_actualiza['fecha_upc_vrs']
        
    return df_filtrado_any, df_filtrado_LRTI, df_filtrado_vrs, df_filtrado_upc

def filtros_IH(df,dias=dias, meses=meses, semana_lower=28,semana_upper=42, con_reemplazo=True):

    #MUERTOS
    if 'VIVO' in(df.columns):
        print("Datos perdidos por muertes: ", df.query('VIVO=="NO"').shape[0])
        df_vivos = df.copy().query('VIVO=="SI"')
        df_muertos = df.copy().query('VIVO=="NO"')
    
    #SEMANA Y PESO GESTACION
    df_vivos_f = df_vivos[(df_vivos['SEMANAS']>=semana_lower) & (df_vivos['SEMANAS']<=semana_upper) &
                      (df_vivos['PESO'] >= df_vivos['p_00001_lognormal']) & (df_vivos['PESO'] <= df_vivos['p_99999_lognormal'])]
    
    # print("Datos perdidos por filtro peso: " , df.drop_duplicates().shape[0]-df_filtered.drop_duplicates().shape[0])
    
    print("Datos perdidos por filtro semanas y peso: " , df_vivos.drop_duplicates().shape[0]-df_vivos_f.drop_duplicates().shape[0])
    
    #drop intersex
    df1 = df_vivos_f.query('SEXO!=9')
    print('Droped intersex:', df_vivos_f.query('SEXO==9').shape[0])
    
    #EDAD MADRE ATIPICA
    print("Datos perdidos por edad madre atípica:", (df1.shape[0] - df1[(df1['EDAD_M']<=51) & (df1['EDAD_M']>=12)].shape[0]))
    df1 = df1[(df1['EDAD_M']<=51) & (df1['EDAD_M']>=12)]
    
    #FECHA ING < FECHA NAC
    df1_outIngNac = df1[(df1['fechaIng_any'] < df1['fecha_nac'])]
    print("Datos perdidos por fecha ingreso menor a fecha nacimiento:", df1_outIngNac.shape[0])
    df1 = df1[~df1.RUN.isin(df1_outIngNac.RUN.unique())]
    
    
    ruts_inmune_precampain = df1.loc[df1.fechaInm < pd.to_datetime("2024-04-01")].RUN.unique()
    df1 =  df1[~df1.RUN.isin(ruts_inmune_precampain)]
    
    #FECHA VRS pre camapaña
    # print("Datos perdidos por vrs pre camapaña:", df1[(df1.vrs_pre_campaña==0)].shape[0])
    # df1 = df1[(df1.vrs_pre_campaña==0)]
    # df1 = df1[~df1.RUN.isin(df1_outIngNac.RUN.unique())]
    
    
    #FECHA INGRESO < FECHA INMUNE
    fechas_precampain = {'fechaIng_vrs': 'vrs_pre_campaña', 'fecha_upc_vrs': 'upc_pre_campaña', 'fechaIng_LRTI': 'lrti_pre_campaña', 'fechaIng_any': 'any_pre_campaña'}
    fechas_df = {'fechaIng_vrs': df1.copy(), 'fecha_upc_vrs': df1.copy(), 'fechaIng_LRTI': df1.copy(), 'fechaIng_any': df1.copy()}
    
    
    
    for key, value in fechas_precampain.items():
        
        fechas_df[key].loc[fechas_df[key][value] == 1, key] = pd.NaT
        
        fechas_df[key].loc[((fechas_df[key][key] - fechas_df[key]['fecha_nac']).dt.days <= 7) , key] = pd.NaT  ## NEW NET ALL CAUSE TESTING ## & (fechas_df[key][value] == 0)
        
        fechas_df[key].loc[(fechas_df[key]['fechaInm'] > fechas_df[key][key]) , 'fechaInm'] = pd.NaT #si es pre_campaña quiero ver cuando se inmuniza & (fechas_df[key][value] == 0)
        
        
        
        if con_reemplazo:
            #if key != 'fechaIng_any':
            print(key, 'Reemplazos n/a net 7 days inmunizado: ', fechas_df[key].loc[((fechas_df[key][key] - fechas_df[key]['fechaInm']).dt.days <= 7) ].shape[0]) #& (fechas_df[key][value] == 0)
            fechas_df[key].loc[((fechas_df[key][key] - fechas_df[key]['fechaInm']).dt.days <= 7)  , 'fechaInm'] = pd.NaT  #si es pre_campaña quiero ver cuando se inmuniza& (fechas_df[key][value] == 0) 
            #else: 
                #fechas_df[key].loc[((fechas_df[key]['fechaIng_LRTI'] - fechas_df[key]['fechaInm']).dt.days <= 7) & (fechas_df[key][value] == 0)  , 'fechaInm'] = pd.NaT
        else:
            fechas_df[key]['between_ingVrs_inm'] = (fechas_df[key][key] - fechas_df[key]['fechaInm']).dt.days
    
    
    #VRS antes del 1 de abril
    df_precampaña = fechas_df['fechaIng_vrs'].copy().query('vrs_pre_campaña==1').sort_values(by='fechaIng_vrs').drop_duplicates(subset=['RUN'], keep='first')
    fechas_df['fechaIng_vrs']['vrs_preabril'] = np.where((fechas_df['fechaIng_vrs'].vrs_pre_campaña==0) & (fechas_df['fechaIng_vrs'].RUN.isin(df_precampaña.RUN.unique())), 1, 0)
    
    #BORRAR DUPLICADOS
    df_aux = fechas_df['fecha_upc_vrs'].copy()
    #df_aux.loc[df_aux.upc_pre_campaña == 1, 'fecha_upc_vrs'] = pd.NaT
    df_upc = df_aux.copy()[df_aux.fecha_upc_vrs.notna()].sort_values(by='fecha_upc_vrs').drop_duplicates(subset=['RUN'], keep='first')
    df_noupc = df_aux.copy()[~df_aux.RUN.isin(df_upc.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_upc = pd.concat([df_noupc,df_upc])
    
    df_aux = fechas_df['fechaIng_vrs'].copy()
    #df_aux.loc[df_aux.vrs_pre_campaña == 1, 'fechaIng_vrs'] = pd.NaT
    vrs = df_aux.copy()[df_aux.fechaIng_vrs.notna()].sort_values(by='fechaIng_vrs').drop_duplicates(subset=['RUN'], keep='first')
    df_novrs = df_aux.copy()[~df_aux.RUN.isin(vrs.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_vrs = pd.concat([df_novrs,vrs])
    
    df_aux = fechas_df['fechaIng_LRTI'].copy()
    #df_aux.loc[df_aux.lrti_pre_campaña == 1, 'fechaIng_LRTI'] = pd.NaT
    LRTI = df_aux.copy()[df_aux.fechaIng_LRTI.notna()].sort_values(by='fechaIng_LRTI').drop_duplicates(subset=['RUN'], keep='first')
    df_noLRTI = df_aux.copy()[~df_aux.RUN.isin(LRTI.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_LRTI = pd.concat([df_noLRTI,LRTI])
    
    df_aux = fechas_df['fechaIng_any'].copy()
    #df_aux.loc[df_aux.DIAG1.str[0] == 'P', 'fechaIng_any'] = pd.NaT
    #df_aux.loc[df_aux.any_pre_campaña == 1, 'fechaIng_any'] = pd.NaT
    df_any = df_aux.copy()[df_aux.fechaIng_any.notna()].sort_values(by='fechaIng_any').drop_duplicates(subset=['RUN'], keep='first')
    df_noany = df_aux.copy()[~df_aux.RUN.isin(df_any.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado_any = pd.concat([df_noany,df_any])
    
    #Crear cumplemeses vivo y cumplesemana inmune
    for df_f in [df_filtrado_vrs, df_filtrado_upc, df_filtrado_LRTI, df_filtrado_any]:
        
        for d in dias:
            nombre_columna = f'inm_{d}d'
            df_f[nombre_columna] = df_f['fechaInm'] + pd.DateOffset(days=d)
            
        for m in meses:
            nombre_columna = f'age_{m}m'
            df_f[nombre_columna] = df_f['fecha_nac'] + pd.DateOffset(months=m)
        
        df_f['inmunizado'] = df_f['fechaInm'].notna().astype(int)
    
    return df_filtrado_any, df_filtrado_LRTI, df_filtrado_vrs, df_filtrado_upc, df_muertos, df_precampaña

def cortes(df_f,dias=dias, meses=meses, cortes = {17:'2024-04-21',18:'2024-04-28',19:'2024-05-05',20:'2024-05-12',
              21:'2024-05-19',22:'2024-05-26',23:'2024-06-02',24:'2024-06-09',
              25:'2024-06-16',26:'2024-06-23',27:'2024-06-30',28:'2024-07-07',
              29:'2024-07-14',30:'2024-07-21',31:'2024-07-28',32:'2024-08-04',
              33:'2024-08-11',34:'2024-08-18',35:'2024-08-25',36:'2024-09-01',
              37:'2024-09-08' ,38:'2024-09-15',39:'2024-09-22',40:'2024-09-29',
              41:'2024-10-6',42:'2024-10-13',43:'2024-10-20'
              },specific_week=None):

    dfs = {}
    
    fechas_events = {'fechaInm':'inmunizado', 'fecha_upc_vrs':'event_upc','fechaIng_vrs':'event_vrs',
                         'fechaIng_vrs_Dall':'event_vrs_Dall', 'fecha_upc_vrs_Dall': 'event_upc_Dall', 
                         'fechaIng_LRTI':'event_LRTI', 'fechaIng_any':'event_any'}
    
    if specific_week:
        df_corte = df_f.copy()
        
        fecha = cortes[specific_week]
        
        fecha_dt = pd.to_datetime(fecha) + pd.DateOffset(days=7)
        
        for m in meses:
            nombre_columna = f'age_{m}m'
            nombre_columna_bin = f'si_{m}_meses'
            df_corte.loc[(df_corte[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
            df_corte[nombre_columna_bin] = 1 - df_corte[nombre_columna].isna()
        
        for d in dias:
            nombre_columna = f'inm_{d}d'
            nombre_columna_bin = f'inm_mayor_{d}d'
            df_corte.loc[(df_corte[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
            df_corte[nombre_columna_bin] = 1 - df_corte[nombre_columna].isna()
        
        
        for fecha_col, event_col in fechas_events.items():
            df_corte.loc[df_corte[fecha_col] > fecha_dt, fecha_col] = pd.NaT
            df_corte[event_col] = df_corte[fecha_col].notnull().astype(int)
        
        df_specific = df_corte.loc[(df_corte.fecha_nac.dt.isocalendar().week <= specific_week) | (df_corte.fecha_nac.dt.isocalendar().year == 2023)]
    
        return df_specific
    
    else:
        for i, fecha in  cortes.items():
            
            df_corte = df_f.copy()
            
            fecha_dt = pd.to_datetime(fecha) + pd.DateOffset(days=7)
            
            for m in meses:
                nombre_columna = f'age_{m}m'
                nombre_columna_bin = f'si_{m}_meses'
                df_corte.loc[(df_corte[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
                df_corte[nombre_columna_bin] = 1 - df_corte[nombre_columna].isna()
            
            for d in dias:
                nombre_columna = f'inm_{d}d'
                nombre_columna_bin = f'inm_mayor_{d}d'
                df_corte.loc[(df_corte[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
                df_corte[nombre_columna_bin] = 1 - df_corte[nombre_columna].isna()
            
            for fecha_col, event_col in fechas_events.items():
                df_corte.loc[df_corte[fecha_col] > fecha_dt, fecha_col] = pd.NaT
                df_corte[event_col] = df_corte[fecha_col].notna().astype(int)
            
            dfs[i] = df_corte.loc[(df_corte.fecha_nac.dt.isocalendar().week <= i) | (df_corte.fecha_nac.dt.isocalendar().year == 2023)]
        
        return dfs

def add_covariate_time(df_0, df_model, fecha_cov, binaria_cov, event, stop_max):
    cv_m = df_model[df_model.fecha == fecha_cov][['RUN','duration',binaria_cov]]
    cv_m = cv_m.rename(columns={'duration':'time'})
    cv_m = cv_m.dropna()
    df_01 = add_covariate_to_timeline(df_0, cv_m, duration_col="time", id_col="RUN", event_col=event) 
    df_01[binaria_cov].fillna(0,inplace=True)
    df_01['inmunizado'].fillna(0,inplace=True)
    df_01['stop'].replace(0,stop_max,inplace=True) #df_0.stop.max() ################################## RRRREEEEEEVISAAAAAR LINEAAAAAAAAAA ######################
    return df_01

def call_data_cox(path,til_week,group_age=False,weeks_inm=False):
    
    dias = list(range(7, (til_week-13)*7, 7)) 
    meses = [1,2,3,4,5,6]

    df_pf = pre_filtred(df_name=path)
    #df_f_any, df_f_LRTI, df_f_vrs, df_f_upc, df_muertos, df_precampaña = filtros_IH(df_pf,dias=dias, meses=meses)
    
    _, _, df_f_vrs, df_f_upc = filtros_IH_new(df_pf,dias=dias, meses=meses)
    
    #df_vrs_tilweek = cortes(df_f_vrs, dias=dias, meses=meses,specific_week=til_week)
    #df_upc_tilweek = cortes(df_f_upc, dias=dias, meses=meses,specific_week=til_week)
    
    df_vrs_tilweek = df_f_vrs.copy()
    df_upc_tilweek = df_f_upc.copy()

    ids_difference = df_vrs_tilweek.set_index('RUN').sort_values(by='RUN').compare(df_upc_tilweek.set_index('RUN').sort_values(by='RUN')).index
    df_f_vrs_diff = df_vrs_tilweek[df_vrs_tilweek.RUN.isin(ids_difference)]
    df_f_upc_diff = df_upc_tilweek[df_upc_tilweek.RUN.isin(ids_difference)]
    df_f_comun = df_vrs_tilweek[~df_vrs_tilweek.RUN.isin(ids_difference)] # Da lo mismos usar vrs o upc, pq son identicas las filas        

    subset = ['RUN',
          "fecha_nac","fechaIng_any","fechaEgr","FECHA_INMUNIZACION","fechaInm","fechaIng_vrs","fechaIng_vrs_Dall", #FECHAS
          "fecha_upc","fecha_upc_vrs","fecha_upc_vrs_Dall", 'fechaIng_LRTI', 'DIAG1', 'MARCA', 'PREVI','DIAG2','DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9','DIAG10', 'DIAG11',
          "VRS_D1","VRS_D1y3","VRS_Dall",'inmunizado', #VRS E INMINE
          "group","sexo","prematuro_extremo","muy_prematuro","prematuro_moderado","prematuro","region","Macrozona2", #covs birth
          'EDAD_M','INS_C_M','REG_RES','URBA_RURAL','COMUNA','COMUNA_N', #covs birth
          'atypic_mom_age', 'INS_N_M', 'SEMANAS', #covs birth #,'vrs_preabril'
          "cama","days_upc","dias_en_ing","days_estad_vrs","days_estad_vrs_Dall", #UPC
          "is_rural","categori_macro","categori_regions","exp_rural","percent_poor","percent_poor_multidim","is_poor","vrs_pre_campaña",'porcent_rural', #otras
          'days_upc_vrs','days_upc_vrs_Dall','event_upc','event_upc_Dall','event_vrs','event_vrs_Dall','event_any','event_LRTI']

    for m in meses:
        nombre_columna = f'age_{m}m'
        nombre_columna_bin = f'si_{m}_meses'
        subset.append(nombre_columna)
        subset.append(nombre_columna_bin)
        
    for d in dias:
        nombre_columna = f'inm_{d}d'
        nombre_columna_bin = f'inm_mayor_{d}d'
        subset.append(nombre_columna)
        subset.append(nombre_columna_bin)
        
    df_vrs_cox = (df_f_vrs_diff.copy()
                  .assign(fecha_evento = lambda x: x.fechaIng_vrs,
                          evento = lambda x: x.event_vrs)
    )
    
    df_upc_cox = (df_f_upc_diff.copy()
                  .assign(fecha_evento = lambda x: x.fecha_upc_vrs,
                          evento = lambda x: x.event_upc)
    )
    
    ################################################################################### PREPROCES EN COMUN #####################################################################

    base_df = df_f_comun.copy()[subset]
    T_inicial = pd.to_datetime('2024-03-31')
    #T_inicial = pd.to_datetime('2023-11-01')
    print('IN preprocess comun')
    
    for m in meses:
        nombre_columna = f'age_{m}m'
        base_df[nombre_columna] = (base_df[nombre_columna] - T_inicial).dt.days
        base_df[nombre_columna] = np.where(base_df[nombre_columna]<=0, 0, base_df[nombre_columna])
        
    for d in dias:
        nombre_columna = f'inm_{d}d'
        base_df[nombre_columna] = (base_df[nombre_columna] - T_inicial).dt.days
        base_df[nombre_columna] = np.where(base_df[nombre_columna]<=0, 0, base_df[nombre_columna])

    fechas_var = ["fecha_nac","fechaInm","fechaIng_vrs","fechaIng_any","fechaIng_LRTI","fechaIng_vrs_Dall","fecha_upc_vrs","fecha_upc_vrs_Dall"]

    for col in fechas_var:
        base_df[col] = (base_df[col] - T_inicial).dt.days

    base_df['start'] = np.where(base_df['fecha_nac']<=0, 0, base_df['fecha_nac'])
    
    for m in meses:
        nombre_columna = f'age_{m}m'
        fechas_var.append(nombre_columna)

    for d in dias:
        nombre_columna = f'inm_{d}d'
        fechas_var.append(nombre_columna)

    df_model = pd.melt(base_df, id_vars=['RUN'], value_vars=fechas_var,
                        var_name='fecha', value_name='duration',ignore_index=True)

    covariates = ['RUN',"VRS_D1","VRS_D1y3","VRS_Dall",'inmunizado',
                        "group","sexo","prematuro_extremo","muy_prematuro","prematuro_moderado","prematuro","region","Macrozona2", 
                        'EDAD_M','INS_C_M','REG_RES','URBA_RURAL','COMUNA','COMUNA_N', 'DIAG1','MARCA', 'PREVI','DIAG2','DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9','DIAG10', 'DIAG11',
                        'atypic_mom_age', 'INS_N_M', 'SEMANAS', #'vrs_preabril',
                        "cama","days_upc","dias_en_ing","days_estad_vrs","days_estad_vrs_Dall", 
                        "is_rural","categori_macro","categori_regions","exp_rural","percent_poor","percent_poor_multidim","is_poor","vrs_pre_campaña",'porcent_rural',
                        'days_upc_vrs','days_upc_vrs_Dall','event_upc','event_upc_Dall','event_vrs','event_vrs_Dall','event_any','event_LRTI']

    for m in meses:
        nombre_columna_bin = f'si_{m}_meses'
        covariates.append(nombre_columna_bin)
        
    for d in dias:
        nombre_columna_bin = f'inm_mayor_{d}d'
        covariates.append(nombre_columna_bin)

    df_model = df_model.merge(base_df[covariates], on = 'RUN' , how='left').assign(group = lambda x: np.where(x['group']=='CATCH_UP',1,0))
    
    subset2 = ['RUN',"group","sexo","prematuro_extremo","muy_prematuro","prematuro_moderado",
           "prematuro","region","Macrozona2",'EDAD_M','INS_C_M','REG_RES','DIAG1','MARCA', 'PREVI','DIAG2','DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9','DIAG10', 'DIAG11',
           'URBA_RURAL','COMUNA','COMUNA_N','atypic_mom_age', 'INS_N_M', 'SEMANAS',"cama",
           "is_rural","categori_macro","categori_regions","exp_rural","percent_poor",#'vrs_preabril',
           "percent_poor_multidim","is_poor","vrs_pre_campaña",'porcent_rural',
           'event_upc','event_upc_Dall','event_vrs', 'event_vrs_Dall','event_any','event_LRTI','duration']

    base_na = df_model[(df_model.fecha=='fechaIng_vrs') & (df_model.duration.isna())][subset2].drop_duplicates(subset=['RUN'])
    base_vrs = df_model[(df_model.fecha=='fechaIng_vrs') & (df_model.duration.notna())][subset2].drop_duplicates(subset=['RUN'])
    base_upc = df_model[(df_model.fecha=='fecha_upc_vrs') & (df_model.duration.notna())][subset2].drop_duplicates(subset=['RUN'])
    
    
    base_na['duration'] = base_na['duration'].fillna(df_model.duration.max())
    base_na = to_long_format(base_na, duration_col="duration")

    base_vrs['duration'] = base_vrs['duration'].fillna(df_model.duration.max())
    base_vrs = to_long_format(base_vrs, duration_col="duration")

    base_upc['duration'] = base_upc['duration'].fillna(df_model.duration.max())
    base_upc = to_long_format(base_upc, duration_col="duration")
    
    subset2.remove("duration")
    subset2_sinevent=[var for var in subset2 if not var.startswith("event_")]
    
    covs_fijas_na = base_na.copy()[subset2_sinevent]
    covs_fijas_vrs = base_vrs.copy()[subset2_sinevent]
    covs_fijas_upc = base_upc.copy()[subset2_sinevent]
    covs_comunes = [covs_fijas_na, covs_fijas_vrs, covs_fijas_upc]
    
    
    base_na = base_na[['RUN','stop','event_vrs']].merge(base_df[['RUN','start']],on='RUN',how='left') 
    base_vrs = base_vrs[['RUN','stop','event_vrs']].merge(base_df[['RUN','start']],on='RUN',how='left') 
    base_upc = base_upc[['RUN','stop','event_upc']].merge(base_df[['RUN','start']],on='RUN',how='left')
    
    stop_max = pd.concat([base_na['stop'], base_vrs['stop'], base_upc['stop']]).max()
    
    df_0_na = add_covariate_time(base_na, df_model, 'fechaInm', 'inmunizado', 'event_vrs',stop_max)
    df_0_vrs = add_covariate_time(base_vrs, df_model, 'fechaInm', 'inmunizado', 'event_vrs',stop_max)
    df_0_upc = add_covariate_time(base_upc, df_model, 'fechaInm', 'inmunizado', 'event_upc',stop_max)
    
    
    if group_age:
        for m in meses:
            fecha_m = f'age_{m}m'
            binaria_m = f'si_{m}_meses'
            df_0_na = add_covariate_time(df_0_na, df_model, fecha_m, binaria_m, 'event_vrs',stop_max)
            df_0_vrs = add_covariate_time(df_0_vrs, df_model, fecha_m, binaria_m, 'event_vrs',stop_max)
            df_0_upc = add_covariate_time(df_0_upc, df_model, fecha_m, binaria_m, 'event_upc',stop_max)
    
    if weeks_inm:
        for d in dias: 
            fecha_d = f'inm_{d}d'
            binaria_d = f'inm_mayor_{d}d'
            df_0_na = add_covariate_time(df_0_na, df_model, fecha_d, binaria_d, 'event_vrs',stop_max)
            df_0_vrs = add_covariate_time(df_0_vrs, df_model, fecha_d, binaria_d, 'event_vrs',stop_max)
            df_0_upc = add_covariate_time(df_0_upc, df_model, fecha_d, binaria_d, 'event_upc',stop_max)
    
    df_0_list = [df_0_na,df_0_vrs,df_0_upc]
    
    ################################################################################### PREPROCES VRS Y UPC #####################################################################
    dfs_cox = [df_vrs_cox, df_upc_cox]
    subset.append('evento')
    subset.append('fecha_evento')
    covariates.append('evento')
    subset2.append('evento')
    
    print('IN preprocess VRS Y UPC')
    for i, df_fe in enumerate(dfs_cox):
        fechas_var = ["fecha_nac","fechaInm","fechaIng_vrs",'fechaIng_any','fechaIng_LRTI',"fechaIng_vrs_Dall","fecha_upc_vrs","fecha_upc_vrs_Dall",'fecha_evento']
        subset2.append('duration')
        base_df = df_fe.copy()[subset]

        for m in meses:
            nombre_columna = f'age_{m}m'
            base_df[nombre_columna] = (base_df[nombre_columna] - T_inicial).dt.days
            base_df[nombre_columna] = np.where(base_df[nombre_columna]<=0, 0, base_df[nombre_columna])
            
        for d in dias:
            nombre_columna = f'inm_{d}d'
            base_df[nombre_columna] = (base_df[nombre_columna] - T_inicial).dt.days
            base_df[nombre_columna] = np.where(base_df[nombre_columna]<=0, 0, base_df[nombre_columna])
             
        for col in fechas_var:
            base_df[col] = (base_df[col] - T_inicial).dt.days

        base_df['start'] = np.where(base_df['fecha_nac']<=0, 0, base_df['fecha_nac'])
        
        for m in meses:
            nombre_columna = f'age_{m}m'
            fechas_var.append(nombre_columna)

        for d in dias:
            nombre_columna = f'inm_{d}d'
            fechas_var.append(nombre_columna)

        df_model = pd.melt(base_df, id_vars=['RUN'], value_vars=fechas_var,
                            var_name='fecha', value_name='duration',ignore_index=True)

        df_model = df_model.merge(base_df[covariates], on = 'RUN' , how='left').assign(group = lambda x: np.where(x['group']=='CATCH_UP',1,0))
        
        base_df_model = df_model[(df_model.fecha=='fecha_evento')][subset2].drop_duplicates(subset=['RUN'])
        
        base_df_model['duration'] = base_df_model['duration'].fillna(df_model.duration.max())
        base_df_model = to_long_format(base_df_model, duration_col="duration")

        subset2.remove("duration")
        subset2_sinevent=[var for var in subset2 if not var.startswith("event")]
        covs_fijas = base_df_model.copy()[subset2_sinevent]

        base_df_model = base_df_model[['RUN','stop','evento']].merge(base_df[['RUN','start']],on='RUN',how='left') 
        
        stop_max = max(base_df_model.stop.max(), stop_max)

        
        df_0 = add_covariate_time(base_df_model, df_model, 'fechaInm', 'inmunizado', 'evento', stop_max)
        
        if group_age:
            for m in meses:
                fecha_m = f'age_{m}m'
                binaria_m = f'si_{m}_meses'
                df_0 = add_covariate_time(df_0, df_model, fecha_m, binaria_m, 'evento', stop_max) #"df_0_na =" edited 04-01-2025
        
        if weeks_inm:
            for d in dias: 
                fecha_d = f'inm_{d}d'
                binaria_d = f'inm_mayor_{d}d'
                df_0 = add_covariate_time(df_0, df_model, fecha_d, binaria_d, 'evento', stop_max) #"df_0_na =" edited 04-01-2025
        
        if i==0:
            df_0 = pd.concat([df_0.rename(columns={'evento': 'event_vrs'}), df_0_list[0], df_0_list[i+1]])
            
        elif i==1:
            df_0 = pd.concat([df_0.rename(columns={'evento': 'event_upc'}), df_0_list[0].copy().rename(columns={'event_vrs': 'event_upc'}), df_0_list[i+1]])
            
        covs = pd.concat([covs_fijas, covs_comunes[0], covs_comunes[1]])
        
        df_0 = (df_0.merge(covs, on='RUN',how='left'))
        
        if i==0:
            df_0_vrs_final = df_0.copy()
            
        elif i==1:
            df_0_upc_final = df_0.copy()

    return df_0_vrs_final, df_0_upc_final, df_vrs_tilweek, df_upc_tilweek # ,df_muertos, df_precampaña

# Esta función, elmina vrs antes de abril y crea variables en base a tiempo
def post_proces_df_cox(df_0):
    if 'start' not in(df_0.columns):
        raise ValueError("Dateframe is not long format")
    
    def assign_month(start_day):
        month_index = int(start_day // 28) % len(months_list)
        return months_list[month_index]

    print('Run nunique stop<start', df_0.query('stop<start').RUN.nunique())
    vrs_preabril = df_0.copy().query('stop<start')
    print('Run nunique start<0', df_0.query('start<0').RUN.nunique())
    
    d = (df_0
         .copy()
         .query('0<=start<=stop')
         .assign(week = lambda x: x.stop // 7, 
                 inmu_time_week = lambda x: x.inmunizado * x.week,
                 months = lambda x: x.stop.apply(assign_month)
                #  group_age = lambda x: x[['si_1_meses', 'si_2_meses', 'si_3_meses',
                #                           'si_4_meses', 'si_5_meses', 'si_6_meses']].apply(lambda row: row.index[row == 1].max(), axis=1).fillna('si_0_meses')
         )
    )

    one_hot_months = pd.get_dummies(d['months'], prefix='month_')

    print("No hay info de probeza de:", d[d.percent_poor_multidim.isna()]['COMUNA'].unique())

    d = (
        d
        .join(one_hot_months)
        .assign(
        mes = lambda x: x['stop'] // 30,
        inmune_time_month = lambda x: x['inmunizado'] * x['mes']
        )
    )
    
    d['group_age'] = d[['si_1_meses', 'si_2_meses', 'si_3_meses', 'si_4_meses', 'si_5_meses', 'si_6_meses']].apply(
        lambda row: row.index[row == 1].max(), axis=1).fillna('si_0_meses')

    d['acumulado_inmune'] = (d.assign(dura_inmu = lambda x: (d['stop'] - d['start']) * d['inmunizado']).groupby('RUN')['dura_inmu'].cumsum())

    d = d.assign(week_being_inmu = lambda x: x.acumulado_inmune // 7,mes_being_inmu = lambda x: x.acumulado_inmune // 28)

    one_hot_weekBI = pd.get_dummies(d['week_being_inmu'], prefix='weekBeingIn')
    d = d.join(one_hot_weekBI)

    one_hot_mesBI = pd.get_dummies(d['mes_being_inmu'], prefix='mesBeingIn')
    d = d.join(one_hot_mesBI)

    one_hot_macr = pd.get_dummies(d['Macrozona2'], prefix='macrozone')
    d = d.join(one_hot_macr)

    one_hot_regions = pd.get_dummies(d['region'], prefix='REGION')
    d = d.join(one_hot_regions)
    
    return d, vrs_preabril, one_hot_weekBI, one_hot_mesBI, one_hot_macr, one_hot_regions

def printSummary(ctv_0):

    sumary = ctv_0.summary
    sumary = (
        sumary
        .assign(effectiveness = lambda x: (1-np.exp(x['coef'])))
        .assign(eff_lower_95 = lambda x: (1-(x['exp(coef) upper 95%'])))
        .assign(eff_upper_95 = lambda x: (1-(x['exp(coef) lower 95%'])))
        .drop(['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','cmp to', 'z','-log2(p)'], axis=1) #,'se(coef)'
        .round(7)
    )
    cols = sumary.columns.tolist()
    cols.insert(1, cols.pop(cols.index('effectiveness')))
    cols.insert(7,cols.pop(cols.index('p')))
    sumary = sumary[cols]
    
    return sumary

def cox(df,covs,prematuros=False):
    df_cox = df.copy()
    if prematuros:
        df_cox = df_cox.query('32<=SEMANAS<=36')
        
    ctv_0 = CoxTimeVaryingFitter()
    strata = []
    if 'region' in covs:
        strata.append('region')
    elif 'Nombre_REGION' in covs:
        strata.append('NOMBRE_REGION')
    if 'group_age' in covs:
        strata.append('group_age')
        
    event = next((cov for cov in covs if cov.startswith('event')), None)
    
    if len(strata)==0:
        ctv_0.fit(df_cox[covs], id_col="RUN", event_col=event, start_col="start", stop_col="stop")
    else:
        ctv_0.fit(df_cox[covs], id_col="RUN", event_col=event, start_col="start", stop_col="stop",strata=strata)
    
    display(printSummary(ctv_0))
    
def cox_return(df,covs,prematuros=False):
    df_cox = df.copy()
    if prematuros:
        df_cox = df_cox.query('SEMANAS<=36') #32<=
        
    ctv_0 = CoxTimeVaryingFitter()
    strata = []
    if 'region' in covs:
        strata.append('region')
    elif 'Nombre_REGION' in covs:
        strata.append('NOMBRE_REGION')
    if 'group_age' in covs:
        strata.append('group_age')
    if 'PREVI' in covs:
        strata.append('PREVI')
    elif 'previ_2' in covs:
        strata.append('previ_2')
    if 'Macrozona2' in covs:
        strata.append('Macrozona2')
    elif 'leo_zonas' in covs:
        strata.append('leo_zonas')
    
        
        
    event = next((cov for cov in covs if cov.startswith('event')), None)
    ################################################################################################################################################################
    if 'sexo' in covs:
        strata.append('sexo')
    ################################################################################################################################################################
    if 'SEMANAS' in covs:
        strata.append('SEMANAS')
    if 'is_rural' in covs:
        strata.append('is_rural')
    if 'prematuro' in covs:
        strata.append('prematuro')
    
    if 'COMUNA' in covs:
        strata.append('COMUNA')
        
    
    if len(strata)==0:
        ctv_0.fit(df_cox[covs], id_col="RUN", event_col=event, start_col="start", stop_col="stop")
    else:
        ctv_0.fit(df_cox[covs], id_col="RUN", event_col=event, start_col="start", stop_col="stop",strata=strata)
    
    return printSummary(ctv_0)

def cox_return_marcel(df,covs,prematuros=False):
    df_cox = df.copy()
    if prematuros:
        df_cox = df_cox.query('SEMANAS<=36') #32<=
        
    ctv_0 = CoxTimeVaryingFitter()
    strata = []
    if 'group_age' in covs:
        strata.append('group_age')
    if 'region' in covs:
        strata.append('region')
    if 'PREVI' in covs:
        strata.append('PREVI')
    elif 'previ_2' in covs:
        strata.append('previ_2')
    if 'Macrozona2' in covs:
        strata.append('Macrozona2')
    elif 'leo_zonas' in covs:
        strata.append('leo_zonas')
    
    event = next((cov for cov in covs if cov.startswith('event')), None)
    
    if len(strata)==0:
        ctv_0.fit(df_cox[covs], id_col="RUN", event_col=event, start_col="start", stop_col="stop")
    else:
        ctv_0.fit(df_cox[covs], id_col="RUN", event_col=event, start_col="start", stop_col="stop",strata=strata)
    
    return printSummary(ctv_0)

def call_data_cox_auxiliar(path,til_week,group_age=False,fecha_vrs='fechaIng_vrs', lrti_name ='LRTI_all_j'):
    
    dias = list(range(7, (til_week-13)*7, 7)) 
    meses = [1,2,3,4,5,6]

    df_pf = pre_filtred(df_name=path, lrti_name = lrti_name)

    #df_f_any, df_f_LRTI, df_f_vrs, df_f_upc, __, __ = filtros_IH(df_pf,dias=dias, meses=meses)
    
    df_f_any, df_f_LRTI, df_f_vrs, df_f_upc = filtros_IH_new(df_pf)
    
    fechas_df = {'fechaIng_vrs': df_f_vrs, 'fecha_upc_vrs': df_f_upc,'fechaIng_LRTI': df_f_LRTI, 'fechaIng_any': df_f_any}
    
    df_f_cox = fechas_df[fecha_vrs].copy()
    
    #df_f_cox_diff = cortes(df_f_cox, dias=dias, meses=meses,specific_week=til_week)
    df_f_cox_diff = df_f_cox.copy()
    
    subset = ['RUN',
          "fecha_nac","fechaIng_any","fechaEgr","FECHA_INMUNIZACION","fechaInm","fechaIng_vrs","fechaIng_vrs_Dall", #FECHAS
          "fecha_upc","fecha_upc_vrs","fecha_upc_vrs_Dall", 'fechaIng_LRTI','DIAG1', 'MARCA','PREVI','DIAG2','DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9','DIAG10', 'DIAG11',
          "VRS_D1","VRS_D1y3","VRS_Dall",'inmunizado', #VRS E INMINE
          "group","sexo","prematuro_extremo","muy_prematuro","prematuro_moderado","prematuro","region","Macrozona2", #covs birth
          'EDAD_M','INS_C_M','REG_RES','URBA_RURAL','COMUNA','COMUNA_N', #covs birth
          'atypic_mom_age', 'INS_N_M', 'SEMANAS', #covs birth
          "cama","days_upc","dias_en_ing","days_estad_vrs","days_estad_vrs_Dall", #UPC
          "is_rural","categori_macro","categori_regions","exp_rural","percent_poor","percent_poor_multidim","is_poor","vrs_pre_campaña",'porcent_rural', #otras
          'days_upc_vrs','days_upc_vrs_Dall','event_upc','event_upc_Dall','event_vrs','event_vrs_Dall','event_any','event_LRTI']

    for m in meses:
        nombre_columna = f'age_{m}m'
        nombre_columna_bin = f'si_{m}_meses'
        subset.append(nombre_columna)
        subset.append(nombre_columna_bin)
        
    for d in dias:
        nombre_columna = f'inm_{d}d'
        nombre_columna_bin = f'inm_mayor_{d}d'
        subset.append(nombre_columna)
        subset.append(nombre_columna_bin)
        
        
    fechas_events = {'fechaInm':'inmunizado', 'fecha_upc_vrs':'event_upc','fechaIng_vrs':'event_vrs',
                    'fechaIng_vrs_Dall':'event_vrs_Dall', 'fecha_upc_vrs_Dall': 'event_upc_Dall', 
                    'fechaIng_LRTI':'event_LRTI', 'fechaIng_any':'event_any'}
    
    event_vrs_var=fechas_events[fecha_vrs]
    
    df_vrs_cox = (df_f_cox_diff.copy()
                  .assign(fecha_evento = lambda x: x[fecha_vrs],
                          evento = lambda x: x[event_vrs_var])
    )

    
    ################################################################################### PREPROCES EN COMUN #####################################################################

    T_inicial = pd.to_datetime('2024-03-31')

    covariates = ['RUN',"VRS_D1","VRS_D1y3","VRS_Dall",'inmunizado','DIAG1', 'DIAG2','DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9','DIAG10', 'DIAG11',
                        "group","sexo","prematuro_extremo","muy_prematuro","prematuro_moderado","prematuro","region","Macrozona2", 
                        'EDAD_M','INS_C_M','REG_RES','URBA_RURAL','COMUNA','COMUNA_N', 
                        'atypic_mom_age', 'INS_N_M', 'SEMANAS', 'MARCA','PREVI',
                        "cama","days_upc","dias_en_ing","days_estad_vrs","days_estad_vrs_Dall", 
                        "is_rural","categori_macro","categori_regions","exp_rural","percent_poor","percent_poor_multidim","is_poor","vrs_pre_campaña",'porcent_rural',
                        'days_upc_vrs','days_upc_vrs_Dall','event_upc','event_upc_Dall','event_vrs','event_vrs_Dall','event_any','event_LRTI']
    
    for m in meses:
        nombre_columna_bin = f'si_{m}_meses'
        covariates.append(nombre_columna_bin)
    
    subset2 = ['RUN',"group","sexo","prematuro_extremo","muy_prematuro","prematuro_moderado",'DIAG1','DIAG2','DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9','DIAG10', 'DIAG11',
           "prematuro","region","Macrozona2",'EDAD_M','INS_C_M','REG_RES','MARCA','PREVI',
           'URBA_RURAL','COMUNA','COMUNA_N','atypic_mom_age', 'INS_N_M', 'SEMANAS',"cama",
           "is_rural","categori_macro","categori_regions","exp_rural","percent_poor",
           "percent_poor_multidim","is_poor","vrs_pre_campaña",'porcent_rural',
           'event_upc','event_upc_Dall','event_vrs', 'event_vrs_Dall','event_any','event_LRTI']

    ################################################################################### PREPROCES VRS Y UPC #####################################################################
    dfs_cox = [df_vrs_cox]
    subset.append('evento')
    subset.append('fecha_evento')
    covariates.append('evento')
    subset2.append('evento')
    
    print('IN preprocess VRS Y UPC')
    for i, df_fe in enumerate(dfs_cox):
        fechas_var = ['fecha_nac', "fechaInm","fechaIng_vrs",'fechaIng_any','fechaIng_LRTI',"fechaIng_vrs_Dall","fecha_upc_vrs","fecha_upc_vrs_Dall",'fecha_evento']
        subset2.append('duration')
        base_df = df_fe.copy()[subset]
        
        for m in meses:
            nombre_columna = f'age_{m}m'
            base_df[nombre_columna] = (base_df[nombre_columna] - T_inicial).dt.days
            base_df[nombre_columna] = np.where(base_df[nombre_columna]<=0, 0, base_df[nombre_columna])
             
        for col in fechas_var:
            base_df[col] = (base_df[col] - T_inicial).dt.days

        base_df['start'] = np.where(base_df['fecha_nac']<=0, 0, base_df['fecha_nac'])
        
        for m in meses:
            nombre_columna = f'age_{m}m'
            fechas_var.append(nombre_columna)

        df_model = pd.melt(base_df, id_vars=['RUN'], value_vars=fechas_var,
                            var_name='fecha', value_name='duration',ignore_index=True)

        df_model = df_model.merge(base_df[covariates], on = 'RUN' , how='left').assign(group = lambda x: np.where(x['group']=='CATCH_UP',1,0))
        
        base_df_model = df_model[(df_model.fecha=='fecha_evento')][subset2].drop_duplicates(subset=['RUN']) #?????drop duplicates?
        
        base_df_model['duration'] = base_df_model['duration'].fillna(df_model.duration.max())
        base_df_model = to_long_format(base_df_model, duration_col="duration")

        subset2.remove("duration")
        subset2_sinevent=[var for var in subset2 if not var.startswith("event")]
        covs_fijas = base_df_model.copy()[subset2_sinevent]

        base_df_model = base_df_model[['RUN','stop','evento']].merge(base_df[['RUN','start']],on='RUN',how='left') 
        
        stop_max = base_df_model.stop.max()

        df_0 = add_covariate_time(base_df_model, df_model, 'fechaInm', 'inmunizado', 'evento', stop_max)
        
        if group_age:
            for m in meses:
                fecha_m = f'age_{m}m'
                binaria_m = f'si_{m}_meses'
                df_0 = add_covariate_time(df_0, df_model, fecha_m, binaria_m, 'evento', stop_max)
        
        df_0_vrs_final = df_0.copy().rename(columns={'evento': event_vrs_var}).merge(covs_fijas, on='RUN',how='left')
        
    return df_0_vrs_final, df_f_cox
    
def filtros_IH_nirsecl(df,dias=dias, meses=meses, semana_lower=24,semana_upper=42):

    if 'VIVO' in(df.columns):
        print("Datos perdidos por muertes: ", df.query('VIVO=="NO"').shape[0])
        df_vivos = df.copy().query('VIVO=="SI"')
        
    n_ruts_ini = df_vivos[(df_vivos['SEMANAS'].notna()) & (df_vivos['PESO'].notna()) & (df_vivos['p_00001_lognormal'].notna()) & (df_vivos['p_99999_lognormal'].notna())].RUN.nunique()
    
    print("ruts perdidos por filtro semanas y peso: " , df_vivos[((df_vivos['SEMANAS']<semana_lower) | (df_vivos['SEMANAS']>semana_upper) |
                                                                   (df_vivos['PESO'] < df_vivos['p_00001_lognormal']) | (df_vivos['PESO'] > df_vivos['p_99999_lognormal']))].RUN.nunique())
    
    #SEMANA Y PESO GESTACION
    df_vivos_f = df_vivos[(df_vivos['SEMANAS']>=semana_lower) & (df_vivos['SEMANAS']<=semana_upper) &
                      (df_vivos['PESO'] >= df_vivos['p_00001_lognormal']) & (df_vivos['PESO'] <= df_vivos['p_99999_lognormal'])]
    
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
    
    ################################################### ANALISIS MARCA ###############################
    
    # print('is_hapenin')
    # media_inmune_catch = np.mean((df1.query('group=="CATCH_UP"').fechaInm - df1.query('group=="CATCH_UP"').fecha_nac).dt.days)
    # media_inmune_nb = np.mean((df1.query('group=="SEASONAL"').fechaInm - df1.query('group=="SEASONAL"').fecha_nac).dt.days)
    
    # df1['random_mark'] = df1['MARCA']

    # marca_1_indices_catchup = df1.query('(group=="CATCH_UP") & (MARCA == 1)').RUN

    # num_to_change = int(0.111 * len(marca_1_indices_catchup))

    # random_run = np.random.choice(marca_1_indices_catchup, size=num_to_change, replace=False)


    # df1.loc[df1.RUN.isin(random_run), 'random_mark'] = 0
    
    # marca_1_indices_nb = df1.query('(group=="SEASONAL") & (MARCA == 1)').RUN

    # num_to_change_nb = int(0.053 * len(marca_1_indices_nb))

    # random_run_nb = np.random.choice(marca_1_indices_nb, size=num_to_change_nb, replace=False)

    # # Cambiar estos índices a 0 en la columna 'random_mark'
    # df1.loc[df1.RUN.isin(random_run_nb), 'random_mark'] = 0
    
    # df1['fechaInm'] = np.where((df1.random_mark==1) & (df1.group=="SEASONAL"), df1.fecha_nac + pd.DateOffset(days = media_inmune_nb), np.where((df1.random_mark==1) & (df1.group=="CATCH_UP"), df1.fecha_nac + pd.DateOffset(days = media_inmune_catch), df1.fechaInm))
    
    ##################################################################################
    
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
    
        for m in meses:
            nombre_columna = f'age_{m}m'
            nombre_columna_bin = f'si_{m}_meses'
            df_f[nombre_columna] = df_f['fecha_nac'] + pd.DateOffset(months=m)
            df_f.loc[(df_f[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
            df_f[nombre_columna_bin] = 1 - df_f[nombre_columna].isna()
        
        for d in dias:
            nombre_columna = f'inm_{d}d'
            nombre_columna_bin = f'inm_mayor_{d}d'
            df_f[nombre_columna] = df_f['fechaInm'] + pd.DateOffset(days=d)
            df_f.loc[(df_f[nombre_columna]> fecha_dt), nombre_columna] = pd.NaT
            df_f[nombre_columna_bin] = 1 - df_f[nombre_columna].isna()
        
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

def pre_filtred_nac_hist(df_name, lrti_name ='LRTI_Flag'):
    
    path_actual = Path.cwd()
    path_data = path_actual.parent / 'Data'
    
    dfpeso = pd.read_csv(path_data / 'gestation_and_weigh_projected.csv', encoding = "latin-1", sep = ",")
    
    df_initial = pd.read_csv(path_data / df_name, encoding = "latin-1", sep = ";")
    print(f'n_rows_inicial= {df_initial.shape[0]}')
    
    #Merging
    df = (df_initial
          .copy()
          .merge(dfpeso[['SEMANAS', 'p_00001_lognormal', 'p_99999_lognormal']], on='SEMANAS',how='left')
    )
    
    cols_diagnostico = ['AREA_FUNC_I','AREAF_1_TRAS', 'AREAF_2_TRAS', 'AREAF_3_TRAS', 'AREAF_4_TRAS', 'AREAF_5_TRAS', 'AREAF_6_TRAS', 'AREAF_7_TRAS', 'AREAF_8_TRAS', 'AREAF_9_TRAS']

    tras_date = {'AREA_FUNC_I': 'fechaIng_any','AREAF_1_TRAS':'fecha_tras_1', 'AREAF_2_TRAS':'fecha_tras_2', 'AREAF_3_TRAS':'fecha_tras_3', 'AREAF_4_TRAS':'fecha_tras_4', 'AREAF_5_TRAS':'fecha_tras_5'
                , 'AREAF_6_TRAS':'fecha_tras_6', 'AREAF_7_TRAS':'fecha_tras_7', 'AREAF_8_TRAS':'fecha_tras_8', 'AREAF_9_TRAS':'fecha_tras_9'}

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
              fechaIng = lambda x: pd.to_datetime({'year': x['ANO_ING'], 'month': x['MES_ING'], 'day': x['DIA_ING']}, format='%Y-%m-%d'),
              fecha_nac = lambda x: pd.to_datetime(x['FECHA_NAC']),
              fechaIng_any = lambda x: pd.to_datetime(x['fechaIng']),
              #fechaEgr = lambda x: pd.to_datetime(x['FECHA_EGRESO']),
            #   FECHA_INMUNIZACION = lambda x: pd.to_datetime(x['FECHA_INMUNIZACION']),
            #   fechaInm = lambda x: pd.to_datetime(x['FECHA_INMUNIZACION']),
              
              
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
              
              
              #features birth baby
              group = lambda x: np.where(x['fecha_nac'] < pd.to_datetime("2024-04-01"), "CATCH_UP", "SEASONAL"),
              sexo = lambda x: np.where(x['SEXO']==1.0,1,0),
              prematuro_extremo = lambda x: np.where((x['SEMANAS']>=22) & (x['SEMANAS']<=27),1,0),
              muy_prematuro = lambda x: np.where((x['SEMANAS']>=28) & (x['SEMANAS']<=32),1,0),
              prematuro_moderado = lambda x: np.where((x['SEMANAS']>=33) & (x['SEMANAS']<=36),1,0),
              prematuro = lambda x: np.where((x['SEMANAS']<=36),1,0),
             # atypic_mom_age = lambda x: np.where((x.EDAD_M>=45) | (x.EDAD_M<=19), 1, 0),
              region = lambda x: x['NOMBRE_REGION'].fillna('DESCONOCIDO').apply(lambda x: mapear_region(x, regiones) if isinstance(x, str) else x),
              Macrozona2 = lambda x: x['region'].map(region_a_macrozona2),
              
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
             # is_rural = lambda x: np.where(x.porcent_rural>0.5,1,0),
              categori_macro = lambda x: x['Macrozona2'].astype('category').cat.codes + 1,
              categori_regions = lambda x: x['region'].astype('category').cat.codes + 1,
             # exp_rural = lambda x: np.exp(x['porcent_rural']),
              #percent_poor = lambda x: x['percent_poor'].fillna(x['percent_poor'].mean()) * 100,
            #  percent_poor_multidim = lambda x: x['percent_poor_multidim'].fillna(x['percent_poor_multidim'].mean()),
             # is_poor = lambda x: np.where(x.percent_poor_multidim >= x.percent_poor_multidim.median(), 1, 0),
              vrs_pre_campaña = lambda x: np.where(x.fechaIng_vrs < pd.to_datetime("2024-04-01"), 1, 0), ################################### '2023-11-01'
              lrti_pre_campaña = lambda x: np.where(x.fechaIng_LRTI < pd.to_datetime("2024-04-01"), 1, 0),
              any_pre_campaña = lambda x: np.where(x.fechaIng_any < pd.to_datetime("2024-04-01"), 1, 0),
              upc_pre_campaña = lambda x: np.where(x.fecha_upc_vrs < pd.to_datetime("2024-04-01"), 1, 0)
              )
    )

    print(np.unique(reg_print))
    
    # DAYS UPC
    # for i in range(1,10):
    #     diff = f'dias_en_area_{i}'
    #     df[diff] = 0
    
    # for index, row in df.copy().query('cama=="UPC"').iterrows():
    #     if pd.isna(row['fecha_tras_1']):
    #         row['days_upc'] = (row['fechaEgr'] - row['fechaIng_any']).days
            
    #     else:
    #         row['dias_en_ing'] = (row['fecha_tras_1'] - row['fechaIng_any']).days
        
    #         for i in range(1, 10):
    #             date_col = f'fecha_tras_{i}'
    #             date_col_next = f'fecha_tras_{i+1}'
    #             diff = f'dias_en_area_{i}'
                
    #             if date_col_next not in tras_date.values() or pd.isna(row[date_col_next]):
    #                 date_col_next = 'fechaEgr'
    #                 row[diff] = (row[date_col_next] - row[date_col]).days
    #                 break
                
    #             row[diff] = (row[date_col_next] - row[date_col]).days

    #         row['days_upc'] += row['dias_en_ing']*row['AREA_FUNC_I']

    #         for i in range(1, 10):
    #             area_col = f'AREAF_{i}_TRAS'
    #             diff = f'dias_en_area_{i}'
    #             if row[area_col]==1:
    #                 row['days_upc'] += row[diff]
            
    #     df.loc[index, 'days_upc'] = np.where(row['days_upc']==0, 1, row['days_upc'])
    #     df.loc[index, 'dias_en_ing'] = row['dias_en_ing']
        
    #     for i in range(1,10):
    #         diff = f'dias_en_area_{i}'
    #         df.loc[index, diff] = row[diff]

    df = df.assign(
       # days_upc_vrs = lambda x: np.where((x['VRS_D1']==1) & (x['days_upc']>0), x['days_upc'], 0),
       # days_upc_vrs_Dall = lambda x: np.where((x['VRS_Dall']==1) & (x['days_upc']>0), x['days_upc'], 0),
        
        #Events
        event_upc = lambda x : (x['fecha_upc_vrs'].notnull()),
        event_upc_Dall = lambda x : (x['fecha_upc_vrs_Dall'].notnull()),
        event_vrs = lambda x : x['fechaIng_vrs'].notnull(),
        event_vrs_Dall = lambda x : x['fechaIng_vrs_Dall'].notnull(),
        event_LRTI = lambda x: x['fechaIng_LRTI'].notnull(),
        event_any = lambda x : x['fechaIng_any'].notnull(),
      #  take_nirse = lambda x : x['fechaInm'].notnull()
        )
    print(f'n_rows_post_prefiltred= {df.shape[0]}')
    
    return df

def filtros_IH_nac(df,dias=dias, meses=meses, semana_lower=24,semana_upper=42):

    # rut_eliminar =  ["bed99009d64eb031ead9235037fc95761d6f334e1d0bc27be4349f1734ca5b2f"]
    # df = df[~df.RUN.isin(rut_eliminar)]
    
    #MUERTOS
    # if 'FECHA_DEF' in df.columns:
    #     df_vivos = df.copy().query('FECHA_DEF.isna()')
    # else: 
    #     df_vivos = df.query('VIVO=="SI"')
    
    df = df.assign(ingreso_mayor_1_mes = lambda x:  (x.fechaIng_any - x.fecha_nac).dt.days)
    
    df['VIVO_mix'] = np.where((df['VIVO']=='SI') | (df['VIVO'].isna()) & (df['FECHA_DEF'].isna()) , 'SI', 'NO')
    
    df_vivos = df.query('(VIVO_mix=="SI") | ((VIVO_mix=="NO") & (ingreso_mayor_1_mes>=7))')
    
    n_ruts_ini = df_vivos[(df_vivos['SEMANAS'].notna()) & (df_vivos['PESO'].notna()) & (df_vivos['p_00001_lognormal'].notna()) & (df_vivos['p_99999_lognormal'].notna())].RUN.nunique()
    
    # print("ruts perdidos por filtro semanas y peso: " , df_vivos[((df_vivos['SEMANAS']<semana_lower) | (df_vivos['SEMANAS']>semana_upper) |
    #                                                                (df_vivos['PESO'] < df_vivos['p_00001_lognormal']) | (df_vivos['PESO'] > df_vivos['p_99999_lognormal']))].RUN.nunique())
    
    #SEMANA Y PESO GESTACION
    df_vivos_f = df_vivos[(df_vivos['SEMANAS']>=semana_lower) & (df_vivos['SEMANAS']<=semana_upper) &
                      (df_vivos['PESO'] > 500) & (df_vivos['PESO'] <= df_vivos['p_99999_lognormal'])] #df_vivos['p_00001_lognormal']
    
    # print("Datos perdidos por filtro peso: " , df.drop_duplicates().shape[0]-df_filtered.drop_duplicates().shape[0])
    
   
    
    #drop intersex
    df1 = df_vivos_f.query('SEXO!=9')
    print('Droped intersex:', df_vivos_f.query('SEXO==9').shape[0])    
    
    
    #FECHA ING < FECHA NAC
    df1_outIngNac = df1[(df1['fechaIng_any'] < df1['fecha_nac'])] # probablemente acá de 0, porque en el filtro de "vivo" ya se tomo la resta de las fechas para filtrar
    print("Datos perdidos por fecha ingreso menor a fecha nacimiento:", df1[df1.RUN.isin(df1_outIngNac.RUN.unique())].RUN.nunique())
    df1 = df1[~df1.RUN.isin(df1_outIngNac.RUN.unique())]
    
    #df1.loc[df1.fechaInm < pd.to_datetime("2024-04-01"), 'fechaInm'] = pd.to_datetime("2024-04-01")
    
    df1 = df1.assign(fechaIng_vrs_copy = lambda x: x.fechaIng_vrs)
    
    df1.loc[df1.fechaIng_any <= df1.fecha_nac + pd.DateOffset(days=7), ['fechaIng_any', 'fechaIng_vrs', 'fecha_upc_vrs','fechaIng_LRTI']] = pd.NaT
    
    # df1 = df1[~df1.RUN.isin(df1[df1.fechaIng_any < pd.to_datetime("2024-04-01")].RUN.unique())]

    print('vrs en los primeros 7 dias de nacer:', df1[df1.RUN.isin(df1[df1['fechaIng_vrs_copy'] <= df1.fecha_nac + pd.DateOffset(days=7)].RUN.unique())].RUN.nunique())
    
    ###############
    df1 = df1[~df1.RUN.isin(df1[df1['fechaIng_vrs_copy'] <= df1.fecha_nac + pd.DateOffset(days=7)].RUN.unique())]
    ###############
    
    n_ruts_fin = df1.RUN.nunique()
    
    print('Ruts eliminados:', n_ruts_ini - n_ruts_fin)
    
    #FECHA INGRESO < FECHA INMUNE
    fechas_df = {'fechaIng_vrs': df1, 'fecha_upc_vrs': df1, 'fechaIng_LRTI': df1, 'fechaIng_any': df1}
    
    # for key in fechas_df.keys():
    #     df_onedit = fechas_df[key].copy()
    #     if key == 'fecha_upc_vrs':
    #         df_onedit = df_onedit.copy()[~df_onedit.RUN.isin(df_onedit[df_onedit['fechaIng_vrs'] < pd.to_datetime("2024-04-01")].RUN.unique())]
    #     else:
    #         df_onedit = df_onedit.copy()[~df_onedit.RUN.isin(df_onedit[df_onedit[key] < pd.to_datetime("2024-04-01")].RUN.unique())]
        
    #     #df_onedit.loc[(df_onedit['fechaInm'] > df_onedit[key]) , 'fechaInm'] = pd.NaT 
    #     #print(key, 'Reemplazos n/a net 7 days inmunizado: ', df_onedit.loc[((df_onedit[key] - df_onedit['fechaInm']).dt.days <= 7)].RUN.nunique()) 
    #     ########################################
    #     #df_onedit.loc[((df_onedit[key] - df_onedit['fechaInm']).dt.days <= 7)  , 'fechaInm'] = pd.NaT  
    #     ########################################
    #     fechas_df[key] = df_onedit
        
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
    
    fechas_events = {'fecha_upc_vrs':'event_upc','fechaIng_vrs':'event_vrs',
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
            df_f[fecha_col] = pd.to_datetime(df_f[fecha_col], infer_datetime_format=True, errors='coerce')
            df_f.loc[df_f[fecha_col] > fecha_dt, fecha_col] = pd.NaT
            df_f[event_col] = df_f[fecha_col].notnull().astype(int)

        df_save = df_f.copy()
        
        # df_save = df_f.loc[(df_f.fecha_nac <= pd.to_datetime("2024-09-30")) | (df_f.fecha_nac.dt.isocalendar().year == 2023)] #2024-09-30

        # # # # df_save = (
        # # # #     df_f
        # # # #     .assign(month_nac = lambda x: x.fecha_nac.dt.month)
        # # # #     .query('month_nac<=9')
        # # # # )
        
        dfs_fecha_actualiza[key] = df_save
        
    df_filtrado_any = dfs_fecha_actualiza['fechaIng_any']
    df_filtrado_LRTI = dfs_fecha_actualiza['fechaIng_LRTI']
    df_filtrado_vrs = dfs_fecha_actualiza['fechaIng_vrs']
    df_filtrado_upc = dfs_fecha_actualiza['fecha_upc_vrs']
        
    return df_filtrado_any, df_filtrado_LRTI, df_filtrado_vrs, df_filtrado_upc

def count_nacs_filter(df,df_filt, filter,pali_ruts=[]):
    
    with open(path_data/'lista_ruts_cardio.pkl', 'rb') as f:
        lista_ruts_cardio = pickle.load(f)

    with open(path_data/'lista_ruts_preterms.pkl', 'rb') as f:
        lista_ruts_preterms = pickle.load(f)
    
    df_nac_per_year = (
        df
        .query(filter)
        .assign(FECHA_NAC = lambda x: pd.to_datetime(x.FECHA_NAC, format = 'mixed'),
                year_nac = lambda x: x.FECHA_NAC.dt.year,
                month_nac = lambda x: x.FECHA_NAC.dt.month,
                elegibilidad_alt = lambda df: (df.year_nac + (df.month_nac >= 10).astype(int)))
        .drop_duplicates('RUN')
        .groupby('elegibilidad_alt',as_index=False)
        .size()
        .rename(columns={'size':'nacs'})
    )

    df_nac_per_year_filt = (
        df_filt
        .query(filter)
        .assign(FECHA_NAC = lambda x: pd.to_datetime(x.FECHA_NAC, format = 'mixed'),
                year_nac = lambda x: x.FECHA_NAC.dt.year,
                month_nac = lambda x: x.FECHA_NAC.dt.month,
                elegibilidad_alt = lambda df: (df.year_nac + (df.month_nac >= 10).astype(int)))
        .drop_duplicates('RUN')
        .groupby('elegibilidad_alt',as_index=False)
        .size()
        .rename(columns={'size':'nacs_filtred'})
        .merge(df_nac_per_year, on='elegibilidad_alt', how='left')
        .assign(lost_filter = lambda x: x.nacs - x.nacs_filtred,
                lost_filter_perc = lambda x: round((x.lost_filter/x.nacs),3))
    )
    return df_nac_per_year_filt



