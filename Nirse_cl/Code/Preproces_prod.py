import pandas as pd
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
warnings.filterwarnings("ignore")

path_actual = Path.cwd()
scope = ["https://spreadsheets.google.com/feeds", "https://www.googleapis.com/auth/spreadsheets",
         "https://www.googleapis.com/auth/drive.file", "https://www.googleapis.com/auth/drive"]

creds = ServiceAccountCredentials.from_json_keyfile_name(path_actual.parent.parent/"credencials.json", scope)
client = gspread.authorize(creds)

spreadsheet_modelos = client.open("Modelos_covariables_nirse")  # Reemplaza con el nombre de tu hoja de cálculo

def son_similares(cadena1, cadena2, umbral=0.6):
    similitud = SequenceMatcher(None, cadena1, cadena2).ratio()
    return similitud >= umbral

def mapear_region(region, regiones_dict, umbral=0.6):
    for key, value in regiones_dict.items():
        if son_similares(region, key, umbral):
            return value
    return None 

dias = [7,14,21,28,35,42,49,56,63]
meses = [1,2,3,4,5,6]


def pre_filtred(df_name,drop_intersex=True,nirse=False):
    
    path_actual = Path.cwd()
    path_data = path_actual.parent / 'Data'
    if nirse:
        path_data = path_actual.parent.parent/'Efectividad_Nirse' / 'Data'
    df = pd.read_csv(path_data / df_name, encoding = "latin-1", sep = ";")
    df['FECHA_NAC'] = pd.to_datetime(df['FECHA_NACIMIENTO'], format='%d%b%Y')
    df['FECHA_INMUNIZACION'] = pd.to_datetime(df['FECHA_INMUNIZACION'], format='%d%b%Y')
    df['FECHA_ING'] = pd.to_datetime(df['FECHA_INGRESO'], format='%d%b%Y')
    df['FECHA_EGR'] = pd.to_datetime(df['FECHA_EGRESO'], format='%d%b%Y')

    diagnosticosVRS = ['J121', 'J205', 'J210','J219', 'B974' ]

    df['VRS'] = np.where(df[['DIAG1'#,'DIAG3','DIAG4','DIAG5','DIAG6','DIAG7','DIAG8','DIAG9','DIAG10','DIAG11'
                            ]].isin(diagnosticosVRS).any(axis=1), 1, 0)

    df['event'] = df['VRS'].astype(bool)
    df['is_inm'] = np.where(((pd.notna(df['FECHA_INMUNIZACION'])) & (df['VRS']==0)) | ((df['FECHA_INMUNIZACION'] + timedelta(days=7) <= df['FECHA_ING']) & (df['VRS']==1)), 1, 0)
    df['group'] = np.where(df['FECHA_NAC']< pd.to_datetime("2024-04-01"), "CATCH_UP", "SEASONAL")
    df['vacunasAlDia'] = np.where(df['VACUNADO']=='SI',1,0)
    
    if drop_intersex==True:
        df = df[df['SEXO']!=9]
        df['SEXO'] = np.where(df['SEXO']==1.0,1,0)

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

    comunas = pd.read_excel(path_data/"comunas.xlsx")
    comunas = comunas.rename(columns = {'C_COM': 'COMUNA_N','NOM_REG':'NOMBRE_REGION'})
    comunas = comunas.drop_duplicates(subset='COMUNA_N')

    diagnosticos_upc = [406, 412, 415, 405, 411, 414]

    cols_diagnostico = ['AREA_FUNC_I','AREAF_1_TRAS', 'AREAF_2_TRAS', 'AREAF_3_TRAS', 'AREAF_4_TRAS', 'AREAF_5_TRAS', 'AREAF_6_TRAS', 'AREAF_7_TRAS', 'AREAF_8_TRAS', 'AREAF_9_TRAS','AREAF_EGR']

    tras_date = {'AREA_FUNC_I': 'FECHA_ING','AREAF_1_TRAS':'fecha_tras_1', 'AREAF_2_TRAS':'fecha_tras_2', 'AREAF_3_TRAS':'fecha_tras_3', 'AREAF_4_TRAS':'fecha_tras_4', 'AREAF_5_TRAS':'fecha_tras_5'
                , 'AREAF_6_TRAS':'fecha_tras_6', 'AREAF_7_TRAS':'fecha_tras_7', 'AREAF_8_TRAS':'fecha_tras_8', 'AREAF_9_TRAS':'fecha_tras_9','AREAF_EGR': 'FECHA_EGR'}

    print(df.shape[0])
    df = df.merge(comunas,how='left',on ='COMUNA_N')
    print("post merge comunas", df.shape[0])
    
    df['NOMBRE_REGION'] = df['NOMBRE_REGION'].fillna('DESCONOCIDO')  
    df['NOMBRE_REGION'] = df['NOMBRE_REGION'].apply(lambda x: mapear_region(x, regiones) if isinstance(x, str) else x)

    df['Macrozona2'] = df['NOMBRE_REGION'].map(region_a_macrozona2)

    for i in range(1, 10):
        year_col = f'ANO_{i}_TRAS'
        month_col = f'MES_{i}_TRAS'
        day_col = f'DIA_{i}_TRAS'
        date_col = f'fecha_tras_{i}'
        df[date_col] = pd.to_datetime({'year': df[year_col], 'month': df[month_col], 'day': df[day_col]}, format='%Y-%m-%d')
        
    for col in cols_diagnostico:
        df[col] = df[col].apply(lambda x: 1 if x in diagnosticos_upc else 0)

    df['cama'] = np.where(df[cols_diagnostico].eq(1).any(axis=1),'UPC', "")

    def obtener_fecha_primer_upc(row):
        for col in cols_diagnostico:
            if row[col] == 1:
                fecha_col = tras_date[col]
                return row[fecha_col]
        return None

    df['fecha_upc'] = df.apply(obtener_fecha_primer_upc, axis=1)
    df['days_upc'] = 0
    df['dias_en_ing'] = (df['fecha_tras_1'] - df['FECHA_ING']).apply(lambda x:  1 + x.days if pd.notna(x) else None)

    for i in range(1, 10):
        date_col = f'fecha_tras_{i}'
        date_col_next = f'fecha_tras_{i+1}'
        diff = f'dias_en_area_{i}'

        if i==9:
            date_col_next = 'FECHA_EGR'
        df[diff] = (df[date_col_next] - df[date_col]).apply(lambda x: 1 + x.days if pd.notna(x) else None)

    df['days_upc'] += df['dias_en_ing'].fillna(0)*df['AREA_FUNC_I']

    for i in range(1, 10):
        area_col = f'AREAF_{i}_TRAS'
        diff = f'dias_en_area_{i}'
        df['days_upc'] += df[diff].fillna(0)*df[area_col]

    df['prematuro_extremo'] = np.where((df['SEMANAS']>=22) & (df['SEMANAS']<=27),1,0)
    df['muy_prematuro'] = np.where((df['SEMANAS']>=28) & (df['SEMANAS']<=32),1,0)
    df['prematuro_moderado'] = np.where((df['SEMANAS']>=33) & (df['SEMANAS']<=36),1,0)
    df['prematuro'] = np.where((df['SEMANAS']<=36),1,0)

    df['days_upc_vrs'] = np.where((df['VRS']==1) & (df['days_upc']>0), df['days_upc'], 0)
    df['days_estad_vrs'] = np.where((df['VRS']==1) & (df['DIAS_ESTAD']>0), df['DIAS_ESTAD'], 0)
    
    df_urba_rural = pd.read_excel(path_data/'urba_rural.xlsx')
    
    poblacion_agrupada = (
        df_urba_rural
        .groupby(['Comuna', 'Area (1=Urbano 2=Rural)'])['Poblacion 2024'].sum().reset_index()
        .rename(columns={'Poblacion 2024': 'Total Poblacion'})
    )


    poblacion_total_comuna = (
        poblacion_agrupada
        .groupby('Comuna')['Total Poblacion'].sum().reset_index()
        .rename(columns={'Total Poblacion': 'Poblacion Total Comuna'})
    )


    df_urba_rural_percent = (
        poblacion_agrupada
        .merge(poblacion_total_comuna, on='Comuna')
        .assign(Porcentaje=lambda x: (x['Total Poblacion'] / x['Poblacion Total Comuna']))
        .pivot(index='Comuna', columns='Area (1=Urbano 2=Rural)', values='Porcentaje')
        .rename(columns={1: 'Porcentaje Urbano', 2: 'Porcentaje Rural'})  # Renombra las columnas correspondientes
        .reset_index()
        .pipe(lambda df: df.rename_axis(None, axis=1))
        .rename(columns={'Comuna':'COMUNA' , 'Porcentaje Rural':'porcent_rural'})
    )

    #df_urba_rural_percent=df_urba_rural_percent.rename(columns={'Comuna':'COMUNA' , 'Porcentaje Rural':'porcent_rural'})
    
    df=df.merge(df_urba_rural_percent, on='COMUNA', how='left')
    print("post merge rural", df.shape[0])

    df['is_rural'] = np.where(df.porcent_rural>0.5,1,0)
    df['atypic_mom_age'] = np.where((df.EDAD_M>=45) | (df.EDAD_M<=19), 1, 0)
    df['categori_macro'] = df['Macrozona2'].astype('category').cat.codes + 1
    df['categori_regions'] = df['NOMBRE_REGION'].astype('category').cat.codes + 1
    df['exp_rural'] = np.exp(df['porcent_rural'])
    

    return df


def filtros_IH(df,dias=dias, meses=meses, semana_lower=28,semana_upper=42, con_reemplazo=True,nirse=False):
    
    path_actual = Path.cwd()
    path_data = path_actual.parent / 'Data'
    if nirse:
        path_data = path_actual.parent.parent/'Efectividad_Nirse' / 'Data'
    dfpeso = pd.read_csv(path_data / 'gestation_and_weigh_projected.csv', encoding = "latin-1", sep = ",")
    df_filtered = df.merge(dfpeso[['SEMANAS', 'p_00001_lognormal', 'p_99999_lognormal']], on='SEMANAS')

    df_filtered = df_filtered[
        (df_filtered['PESO'] >= df_filtered['p_00001_lognormal']) &
        (df_filtered['PESO'] <= df_filtered['p_99999_lognormal'])
    ]
    print("Datos perdidos por filtro peso: " , df.drop_duplicates().shape[0]-df_filtered.drop_duplicates().shape[0])
    df1 = df_filtered[(df_filtered['SEMANAS']>=semana_lower) & (df_filtered['SEMANAS']<=semana_upper)]
    print("Datos perdidos por filtro semanas y peso: " , df.drop_duplicates().shape[0]-df1.drop_duplicates().shape[0])

    df1.loc[(df1['FECHA_INMUNIZACION'] > df1['FECHA_ING']) & (df1.event ==1 ), 'FECHA_INMUNIZACION'] = pd.NaT    
    df1.loc[ (~df1['event']) & (~df1['FECHA_ING'].isna()), 'FECHA_ING'] = pd.NaT
    
    if con_reemplazo:
        df1.loc[((df1['FECHA_ING'] - df1['FECHA_INMUNIZACION']).dt.days <= 7) & (df1['event']) & (~df1.FECHA_INMUNIZACION.isna()), 'FECHA_INMUNIZACION'] = pd.NaT
    else:
        df1['vrs_preinmune'] = np.where((df1.event == 1) & (df1.FECHA_ING < df1.FECHA_INMUNIZACION)) 
        
    df1_outIngNac = df1[(df1['FECHA_ING'] < df1['FECHA_NAC']) & (~df1['FECHA_ING'].isna())]
    print("Datos perdidos por fecha ingreso menor a fecha nacimiento:", df1_outIngNac.shape[0])
    df_2=df1[~df1.RUN.isin(df1_outIngNac.RUN)]
    
    df_filtrado=df_2[(df_2['EDAD_M']<=51) & (df_2['EDAD_M']>=12)]
    print("Datos perdidos por edad madre atípica:", (df_2.shape[0] - df_filtrado.shape[0]))
    
    df_filtrado['inmunizado'] = df_filtrado['FECHA_INMUNIZACION'].notna().astype(int)

    vrs = df_filtrado[df_filtrado.event == 1].drop_duplicates(subset=['RUN'], keep='first')
    df_novrs = df_filtrado[~df_filtrado.RUN.isin(vrs.RUN.unique())].drop_duplicates(subset=['RUN'], keep='first')
    df_filtrado = pd.concat([df_novrs,vrs])
    
    for d in dias:
        nombre_columna = f'inm_{d}d'
        df_filtrado[nombre_columna] = df_filtrado['FECHA_INMUNIZACION'] + pd.DateOffset(days=d)
        
    for m in meses:
        nombre_columna = f'age_{m}m'
        df_filtrado[nombre_columna] = df_filtrado['FECHA_NAC'] + pd.DateOffset(months=m)
    #for bucle above replace:
    #df_filtrado['inm_7d'] = df_filtrado['FECHA_INMUNIZACION'] + pd.DateOffset(days=7)
    #df_filtrado['inm_14d'] = df_filtrado['FECHA_INMUNIZACION'] + pd.DateOffset(days=14)
    #df_filtrado['inm_21d'] = df_filtrado['FECHA_INMUNIZACION'] + pd.DateOffset(days=21)
    #df_filtrado['inm_28d'] = df_filtrado['FECHA_INMUNIZACION'] + pd.DateOffset(days=28)

    #df_filtrado['age_3m'] = df_filtrado['FECHA_NAC'] + pd.DateOffset(months=3)
    #df_filtrado['age_6m'] = df_filtrado['FECHA_NAC'] + pd.DateOffset(months=6)
    
    df_filtrado['sex*prem'] = df_filtrado['muy_prematuro']*df_filtrado['SEXO']
    
    return df_filtrado


def cortes(df_f,dias=dias, meses=meses, cortes = {17:'2024-04-21',18:'2024-04-28',19:'2024-05-05',20:'2024-05-12',
              21:'2024-05-19',22:'2024-05-26',23:'2024-06-02',24:'2024-06-09',
              25:'2024-06-16',26:'2024-06-23',27:'2024-06-30',28:'2024-07-07',
              29:'2024-07-14',30:'2024-07-21',31:'2024-07-28',32:'2024-08-04',
              33:'2024-08-11',34:'2024-08-18',35:'2024-08-25',36:'2024-09-01',
              37:'2024-09-08' ,38:'2024-09-15'#,39:'2024-08-25',40:'2024-08-25',41:'2024-08-25'
              }):
    
    #df_f['id'] = df_f['RUN'] + df_f['FECHA_ING'].astype('str')

    dfs = {}
    
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
            
        df_corte.loc[(df_corte['FECHA_INMUNIZACION']> fecha_dt), 'FECHA_INMUNIZACION'] = pd.NaT
        df_corte['inmunizado'] = df_corte['FECHA_INMUNIZACION'].notna().astype(int)
        
        df_corte.loc[(df_corte['FECHA_ING']> fecha_dt), 'FECHA_ING'] = pd.NaT
        df_corte['event'] = df_corte['FECHA_ING'].notna()
        
        #df_corte = df_corte[((df_corte['FECHA_ING'].isna()) | (df_corte['FECHA_ING']<=fecha_dt)) & ((df_corte['FECHA_INMUNIZACION'].isna()) | (df_corte['FECHA_INMUNIZACION']<=fecha_dt))]
        
        dfs[i] = df_corte.loc[(df_corte.FECHA_NAC.dt.isocalendar().week <= i) | (df_corte.FECHA_NAC.dt.isocalendar().year == 2023)]
    
    return dfs

def preprocess_model_IH(df_f):
    
    #['RUN','SEXO','sex*prem','FECHA_NAC','group',
    #                'prematuro_extremo','muy_prematuro','prematuro_moderado',
    #                'prematuro','FECHA_INMUNIZACION','FECHA_ING','inmunizado',
    #                'event','age_3m','age_6m','inm_7d','inm_14d','inm_21d','inm_28d',
    #                'si_3_meses','si_6_meses','inm_mayor_7d','inm_mayor_14d','inm_mayor_21d','inm_mayor_28d',
    #                'EDAD_M','INS_C_M','COMUNA','COMUNA_N','REG_RES','URBA_RURAL', 'Macrozona2','NOMBRE_REGION']
    
    subset=['RUN','SEXO','sex*prem','FECHA_NAC','group',
            'prematuro_extremo','muy_prematuro','prematuro_moderado','porcent_rural',
            'prematuro','FECHA_INMUNIZACION','FECHA_ING','inmunizado','is_rural',
            'event','EDAD_M','INS_C_M','COMUNA','COMUNA_N','REG_RES','URBA_RURAL', 'Macrozona2','NOMBRE_REGION',
            'atypic_mom_age'
            ]
    
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
    
    base_df = df_f[subset]
    
    base_df = base_df.sort_values(by='event')
    #vrs = base_df[base_df.event == 1]
    #df_vrs = base_df[base_df.RUN.isin(vrs.RUN.unique())]
    #df_vrs = df_vrs.drop_duplicates(subset=['RUN'], keep='first')
    #base_df = base_df[~base_df.RUN.isin(vrs.RUN.unique())]
    #base_df = pd.concat([base_df,df_vrs])
    T_inicial = pd.to_datetime('2024-03-31')
    
    
    for m in meses:
        nombre_columna = f'age_{m}m'
        base_df[nombre_columna] = (base_df[nombre_columna] - T_inicial).dt.days
        base_df[nombre_columna] = np.where(base_df[nombre_columna]<=0, 0, base_df[nombre_columna])
        
    for d in dias:
        nombre_columna = f'inm_{d}d'
        base_df[nombre_columna] = (base_df[nombre_columna] - T_inicial).dt.days
        base_df[nombre_columna] = np.where(base_df[nombre_columna]<=0, 0, base_df[nombre_columna])
    
    
    #base_df['age_3m'] = (base_df['age_3m'] - T_inicial).dt.days
    #base_df['age_6m'] = (base_df['age_6m'] - T_inicial).dt.days
    #base_df['inm_7d'] = (base_df['inm_7d'] - T_inicial).dt.days
    #base_df['inm_14d'] = (base_df['inm_14d'] - T_inicial).dt.days
    #base_df['inm_21d'] = (base_df['inm_21d'] - T_inicial).dt.days
    #base_df['inm_28d'] = (base_df['inm_28d'] - T_inicial).dt.days
    
    base_df['FECHA_INMUNIZACION'] = (base_df['FECHA_INMUNIZACION'] - T_inicial).dt.days
    base_df['FECHA_ING'] = (base_df['FECHA_ING'] - T_inicial).dt.days
    base_df['FECHA_NAC'] = (base_df['FECHA_NAC'] - T_inicial).dt.days
    
    #base_df['start'] = np.where(base_df['group']=="CATCH_UP", 0, base_df['FECHA_NAC'])
    
    base_df['start'] = np.where(base_df['FECHA_NAC']<=0, 0, base_df['FECHA_NAC'])
    
    #base_df['age_3m'] = np.where(base_df['age_3m']<=0, 0, base_df['age_3m'])
    #base_df['age_6m'] = np.where(base_df['age_6m']<=0, 0, base_df['age_6m'])
    #base_df['inm_7d'] = np.where(base_df['inm_7d']<=0, 0, base_df['inm_7d'])
    #base_df['inm_14d'] = np.where(base_df['inm_14d']<=0, 0, base_df['inm_14d'])
    #base_df['inm_21d'] = np.where(base_df['inm_21d']<=0, 0, base_df['inm_21d'])
    #base_df['inm_28d'] = np.where(base_df['inm_28d']<=0, 0, base_df['inm_28d'])
    
    cumpleaños_covs=['FECHA_INMUNIZACION', 'FECHA_ING', 'FECHA_NAC']
    
    for m in meses:
        nombre_columna = f'age_{m}m'
        cumpleaños_covs.append(nombre_columna)

    for d in dias:
        nombre_columna = f'inm_{d}d'
        cumpleaños_covs.append(nombre_columna)
    
    #value_vars=['FECHA_INMUNIZACION', 'FECHA_ING','age_3m','age_6m', 'inm_7d','inm_14d','inm_21d','inm_28d']
    
    df_model = pd.melt(base_df, id_vars=['RUN'], value_vars=cumpleaños_covs,
                       var_name='fecha', value_name='duration',ignore_index=True)
    
    #covariates=['RUN','SEXO','sex*prem','FECHA_NAC','group','prematuro_extremo',
    #            'muy_prematuro','prematuro_moderado','prematuro','inmunizado','event','si_3_meses','si_6_meses',
    #            'inm_mayor_7d','inm_mayor_14d','inm_mayor_21d','inm_mayor_28d','EDAD_M','INS_C_M','COMUNA',
    #            'COMUNA_N','REG_RES','URBA_RURAL', 'Macrozona2','NOMBRE_REGION'
    #            ]
    covariates_fijas = ['RUN','SEXO','sex*prem','group','prematuro_extremo', #'FECHA_NAC'
                  'muy_prematuro','prematuro_moderado','prematuro','inmunizado','event',
                  'EDAD_M','INS_C_M','COMUNA','COMUNA_N','REG_RES','URBA_RURAL', 'Macrozona2',
                  'NOMBRE_REGION','is_rural','porcent_rural','atypic_mom_age'
                  ]
    
    covariates = covariates_fijas.copy()
    
    for m in meses:
        nombre_columna_bin = f'si_{m}_meses'
        covariates.append(nombre_columna_bin)
        
    for d in dias:
        nombre_columna_bin = f'inm_mayor_{d}d'
        covariates.append(nombre_columna_bin)
    
    df_model = df_model.merge(base_df[covariates], on = 'RUN' , how='left')
    
    df_model['group']=np.where(df_model['group']=='CATCH_UP',1,0)

    base = df_model[df_model.fecha=='FECHA_ING'][['RUN','SEXO','sex*prem','group',
                                                  'prematuro_extremo','muy_prematuro','prematuro_moderado',
                                                  'prematuro','duration','event','EDAD_M','INS_C_M','COMUNA',
                                                  'COMUNA_N','REG_RES','URBA_RURAL', 'Macrozona2','NOMBRE_REGION',
                                                  'is_rural','porcent_rural','atypic_mom_age']].drop_duplicates(subset=['RUN'])
    
    base['duration'] = base['duration'].fillna(df_model.duration.max())

    base = to_long_format(base, duration_col="duration")

    base_fp=base
    #covs_fijas1 = base_fp[['RUN','prematuro']]
    #covs_fijas2 = base_fp[['RUN','prematuro_extremo','muy_prematuro','prematuro_moderado']]
    covs_fijas3 = base_fp[['RUN','SEXO','muy_prematuro']]
    #covs_fijas4 = base_fp[['RUN','SEXO','prematuro_extremo','muy_prematuro','prematuro_moderado']]
    #covs_fijas5 = base_fp[['RUN','SEXO','muy_prematuro']]
    covs_fijas_news = base_fp[['RUN','SEXO','group','muy_prematuro','EDAD_M',
                               'INS_C_M','COMUNA','COMUNA_N','REG_RES','URBA_RURAL',
                               'Macrozona2','NOMBRE_REGION','is_rural','porcent_rural','atypic_mom_age']]
    covs = [covs_fijas_news, covs_fijas3]

    base = base[['RUN','stop','event']].merge(base_df[['RUN','start']],on='RUN',how='left') 

    cv = df_model[df_model.fecha=='FECHA_INMUNIZACION'][['RUN','duration','inmunizado']]
    cv = cv.rename(columns={'duration':'time'})
    cv = cv.dropna()
    
    df_0 = add_covariate_to_timeline(base, cv, duration_col="time", id_col="RUN", event_col="event") 

    df_0['inmunizado'].fillna(0,inplace=True)
    df_0['stop'].replace(0,base.stop.max(),inplace=True) ################ base.stop.max() into df_model.duration.max()????
    
    # d['inmu_time'] = d['inmunizado']*(d['stop']-d['start'])
    # bins = list(np.arange(0, 64, 7))  + [float('inf')]
    # labels = [f"{i}-{i+7}" for i in range(0, 63, 7)] + ['63+'] 
    # d['inmu_time_groups'] = pd.cut(d['inmu_time'], bins=bins, labels=labels, right=False, ordered=True)
    # d['inmu_time_strata'] = d['inmu_time_strata'].astype('category').cat.codes + 1
    
    return df_0, base, df_model, covs

def cov_3meses(df_0, base, df_model):
    cv_3m = df_model[df_model.fecha=='age_3m'][['RUN','duration','si_3_meses']]
    cv_3m = cv_3m.rename(columns={'duration':'time'})
    cv_3m = cv_3m.dropna()
    df_01 = add_covariate_to_timeline(df_0, cv_3m, duration_col="time", id_col="RUN", event_col="event") 
    df_01['si_3_meses'].fillna(0,inplace=True)
    df_01['inmunizado'].fillna(0,inplace=True)
    df_01['stop'].replace(0,base.stop.max(),inplace=True)
    return df_01

def cov_Mmeses(df_0, base, df_model, m):
    fecha_m = f'age_{m}m'
    col=f'si_{m}_meses'
    cv_m = df_model[df_model.fecha == fecha_m][['RUN','duration',col]]
    cv_m = cv_m.rename(columns={'duration':'time'})
    cv_m = cv_m.dropna()
    df_01 = add_covariate_to_timeline(df_0, cv_m, duration_col="time", id_col="RUN", event_col="event") 
    df_01[col].fillna(0,inplace=True)
    df_01['inmunizado'].fillna(0,inplace=True)
    df_01['stop'].replace(0,base.stop.max(),inplace=True)
    return df_01

def cov_inmune_7(df_0, base, df_model):
    cv_7in = df_model[df_model.fecha=='inm_7d'][['RUN','duration','inm_mayor_7d']]
    cv_7in = cv_7in.rename(columns={'duration':'time'})
    cv_7in = cv_7in.dropna()
    df_01 = add_covariate_to_timeline(df_0, cv_7in, duration_col="time", id_col="RUN", event_col="event") 
    df_01['inm_mayor_7d'].fillna(0,inplace=True)
    df_01['inmunizado'].fillna(0,inplace=True)
    df_01['stop'].replace(0,base.stop.max(),inplace=True)
    return df_01

def cov_inmune_d(df_0, base, df_model, d):
    cv_in = df_model[df_model.fecha == f'inm_{d}d'][['RUN', 'duration', f'inm_mayor_{d}d']]
    cv_in = cv_in.rename(columns={'duration': 'time'})
    cv_in = cv_in.dropna()
    df_01 = add_covariate_to_timeline(df_0, cv_in, duration_col="time", id_col="RUN", event_col="event") 
    df_01[f'inm_mayor_{d}d'].fillna(0, inplace=True)
    df_01['inmunizado'].fillna(0, inplace=True)
    df_01['stop'].replace(0, base.stop.max(), inplace=True)
    return df_01

def cox_model(df):
    ctv_0 = CoxTimeVaryingFitter()
    ctv_0.fit(df, id_col="RUN", event_col="event", start_col="start", stop_col="stop", show_progress=True)
    ctv_0.print_summary()
    ctv_0.plot()

    coef_0 = ctv_0.params_
    conf_0 = ctv_0.confidence_intervals_

    hazard_ratios_0 = 1-np.exp(coef_0)
    hazard_ratios_conf_int_0 = 1-np.exp(conf_0)
    print(hazard_ratios_0)
    print(hazard_ratios_conf_int_0)
    
    
    
def preprocess_model_upc(df_f):
    
    df_f['event'] = (df_f['fecha_upc'].notna()) & (df_f['VRS']==1)
    
    df_f.loc[ (~df_f['event']), 'fecha_upc'] = pd.NaT
    df_f.loc[(df_f['FECHA_INMUNIZACION'] >= df_f['fecha_upc']) & (df_f.event ==1 ), 'FECHA_INMUNIZACION'] = pd.NaT
    base_df = df_f[['RUN','SEXO','sex*prem','prematuro_extremo',
                    'muy_prematuro','prematuro_moderado','prematuro',
                    'FECHA_NAC','FECHA_INMUNIZACION','FECHA_ING','fecha_upc',
                    'inmunizado','group','event','age_3m','age_6m','si_3_meses','si_6_meses']]
    base_df = base_df.sort_values(by='event')

    vrs = base_df[base_df.event == 1]
    df_vrs = base_df[base_df.RUN.isin(vrs.RUN.unique())]
    df_vrs = df_vrs.drop_duplicates(subset=['RUN'], keep='last')
    base_df = base_df[~base_df.RUN.isin(vrs.RUN.unique())]
    base_df = pd.concat([base_df,df_vrs])

    T_inicial = pd.to_datetime('2024-03-31')

    base_df['age_3m'] = (base_df['age_3m'] - T_inicial).dt.days
    base_df['age_6m'] = (base_df['age_6m'] - T_inicial).dt.days
    base_df['FECHA_INMUNIZACION'] = (base_df['FECHA_INMUNIZACION'] - T_inicial).dt.days
    base_df['FECHA_ING'] = (base_df['FECHA_ING'] - T_inicial).dt.days
    base_df['FECHA_NAC'] = (base_df['FECHA_NAC'] - T_inicial).dt.days
    base_df['fecha_upc'] = (base_df['fecha_upc'] - T_inicial).dt.days

    base_df['start'] = np.where(base_df['group']=="CATCH_UP", 0, base_df['FECHA_NAC'])

    base_df['age_3m'] = np.where(base_df['age_3m']<=0, 0, base_df['age_3m'])
    base_df['age_6m'] = np.where(base_df['age_6m']<=0, 0, base_df['age_6m'])

    df_model = pd.melt(base_df, id_vars=['RUN'], value_vars=['FECHA_INMUNIZACION', 'FECHA_ING','fecha_upc','age_3m','age_6m'],var_name='fecha', value_name='duration',ignore_index=True)
    df_model = df_model.merge(base_df[['RUN','inmunizado','event','si_3_meses','si_6_meses','SEXO','sex*prem','prematuro_extremo','muy_prematuro','prematuro_moderado','prematuro'
                                    ]], on = 'RUN' , how='left')

    base = df_model[df_model.fecha=='fecha_upc'][['RUN','SEXO','sex*prem','prematuro_extremo','muy_prematuro','prematuro_moderado','prematuro','duration','event']].drop_duplicates(subset=['RUN'])
    base['duration'] = base['duration'].fillna(df_model.duration.max()) ############################################ SE PONE LA MAXIMA FECHA POSIBLE

    base = to_long_format(base, duration_col="duration")

    base_fp=base
    covs_fijas1 = base_fp[['RUN','prematuro']]
    covs_fijas2 = base_fp[['RUN','muy_prematuro','prematuro_moderado']]
    covs_fijas3 = base_fp[['RUN','SEXO','muy_prematuro']]
    covs_fijas4 = base_fp[['RUN','SEXO','muy_prematuro','prematuro_moderado']]
    covs_fijas5 = base_fp[['RUN','SEXO','muy_prematuro',"sex*prem"]]
    covs=[covs_fijas1, covs_fijas2, covs_fijas3, covs_fijas4, covs_fijas5]

    base = base[['RUN','stop','event']].merge(base_df[['RUN','start']],on='RUN',how='left')
    cv = df_model[df_model.fecha=='FECHA_INMUNIZACION'][['RUN','duration','inmunizado']]
    cv = cv.rename(columns={'duration':'time'})
    cv = cv.dropna()

    df_upc = add_covariate_to_timeline(base, cv, duration_col="time", id_col="RUN", event_col="event")
    df_upc['inmunizado'].fillna(0,inplace=True)
    df_upc['stop'].replace(0,base.stop.max(),inplace=True)

    return df_upc, base, df_model, covs


def split_df(df, split_col):
    uniques = df[split_col].unique()  
    print(f'Hay {df[split_col].nunique()} valores únicos en {split_col}')
    print(f'Hay {df[split_col].isna().sum()} N/A en {split_col}')
    
    ruta_carpeta_base = Path('../Data') / split_col
    ruta_carpeta_base.mkdir(parents=True, exist_ok=True)

    for i in uniques:
        if pd.isnull(i):
            print("Hay un valor N/A como unique, saltando este valor.")
            continue
        else:

            var_name = f'df_UW_{i}'.replace(" ", "_")
            name_split = f'{var_name}.csv'

            splited = df[df[split_col] == i]
            
            ruta_archivo = ruta_carpeta_base / name_split
            splited.to_csv(ruta_archivo, index=False)
            
            globals()[var_name] = splited

            print(f'Variable {var_name} creada con {len(splited)} filas y guardada en {ruta_archivo}')
            
def printSummary(ctv_0,spreadsheet=spreadsheet_modelos):

    sumary = ctv_0.summary
    sumary = (
        sumary
        .assign(effectiveness = lambda x: (1-np.exp(x['coef'])))
        .assign(eff_lower_95 = lambda x: (1-(x['exp(coef) upper 95%'])))
        .assign(eff_upper_95 = lambda x: (1-(x['exp(coef) lower 95%'])))
        .drop(['exp(coef)','se(coef)','exp(coef) lower 95%','exp(coef) upper 95%','cmp to', 'z','-log2(p)'], axis=1)
        .round(3)
    )
    cols = sumary.columns.tolist()
    cols.insert(1, cols.pop(cols.index('effectiveness')))
    cols.insert(7,cols.pop(cols.index('p')))
    sumary = sumary[cols]
    
    # spreadsheet.add_worksheet(title=nombre, rows="100", cols="20")

    # worksheet = spreadsheet.worksheet("Modelo")

    # model = printSummary(ctv_0)

    # df_values = model.reset_index().values.tolist()

    # # Añadir los encabezados (opcional)
    # df_headers = ['Covariates'] + model.columns.tolist()

    # worksheet.clear()

    # worksheet.update('A1', [df_headers] + df_values)
    
    return sumary