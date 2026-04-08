import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import sys
import os
from pathlib import Path
from dateutil.relativedelta import relativedelta
parent_dir = os.getcwd()
#if parent_dir not in sys.path: sys.path.insert(0,parent_dir+'/src')
# import config as cf #from config 

# folder_path_input = cf.dataDirI2
# folder_path = cf.dataDirO2





####### ICD10 ####### 
def diag_ICD10_cat(d):
    if isinstance(d, str) and d not in [None, 'nan', 'NaN']:
        if d.startswith('J'):
            # Acute Bronchitis/Bronchiolitis J20-J21
            if d.startswith(('J20', 'J21')):
                return 'acute_bronchitis_bronchiolitis'
            # Pneumonia J12-J18
            elif d[1:3] in [f'{x}' for x in range(12,19)]:
                return 'pneumonia'
            # Influenza J09-J11
            elif d.startswith(('J09', 'J10', 'J11')):
                return 'influenza'
            # Upper respiratory tract infection J00-J06
            elif d[1:3] in [f'0{x}' for x in range(0,6)]:
                return 'upper_respiratory_tract_infection'
            # Obtructive Pulmonary Disease J40-J46
            elif d[1:3] in [f'{x}' for x in range(40,47)]:
                return 'obtructive_pulmonary_disease'
            # Other respiratory conditions J22, J30-J39, J47, J60-J98
            elif d.startswith(('J22','J47')) or d[1:3] in [f'{x}' for x in range(30,40)] or d[1:3] in [f'{x}' for x in range(60,99)]:
                return 'other_respiratory_conditions'
            else:
                return 'non_respiratory_disease'
        elif d.startswith(('U071','U072')):
            return 'covid'
    
    return 'non_respiratory_disease'


####### catchup e inseason ####### 
def catchup(fechaNac, y):
    fechas = [pd.Timestamp(year=y, month=m, day=1) for m in range(4, 0, -1)] + \
             [pd.Timestamp(year=y-1, month=m, day=1) for m in range(12, 3, -1)]

    for i, fecha_limite in enumerate(fechas):
        if fechaNac >= fecha_limite:
            return i
    return -1

def catchup2(month_born,year_born,year_in):
    return ((year_in-year_born)==1)*(month_born>=10)*1 + (((year_in-year_born)==0)*(month_born<=3))*1

def inseason(month_born,year_born,year_in):
    return (((year_in-year_born)==0)*(month_born>=4))*(month_born<=9)*1

def inseason2(month_born,year_born,year_in):
    return (((year_in-year_born)==0)&(month_born>=4))&(month_born<=9)

####### Transformar tabla egresos ####### 
def transform_translations(dataframe, areasUPC, areasMB, min_fecha):
    tras_cols = [col for col in dataframe.columns if '_TRAS' in col]
    base_cols = [col for col in dataframe.columns if '_TRAS' not in col]
    steps = sorted({int(col.split('_')[1]) for col in tras_cols})

    areasUPC_set = set(areasUPC)  
    areasMB_set = set(areasMB)  
    rows = []
    l = 0
    for row in dataframe.itertuples(index=False):
        
        new_row = {col: getattr(row, col) for col in base_cols}
        new_row.update({
            'DIA_ING': getattr(row, 'DIA_ING', None),
            'MES_ING': getattr(row, 'MES_ING', None),
            'AÑO_ING': getattr(row, 'ANO_ING', None),
            'AREA_FUNC_I': getattr(row, 'AREA_FUNC_I', None),
            'UPC': int(getattr(row, 'AREA_FUNC_I', None) in areasUPC_set),
            'MB': int(getattr(row, 'AREA_FUNC_I', None) in areasMB_set),
            'SER_CLIN_I': getattr(row, 'SER_CLIN_I', None),
            'fechaIng': getattr(row, 'fechaIng', None)
        })
              
        new_row.update({'indexIng': (new_row.get('fechaIng', None)-min_fecha).days})
        new_row.update({'indice': l})

        last_egr = {
            'DIA_EGR': getattr(row, 'DIA_EGR', None),
            'MES_EGR': getattr(row, 'MES_EGR', None),
            'AÑO_EGR': getattr(row, 'ANO_EGR', None),
            'AREAF_EGR': getattr(row, 'AREAF_EGR', None),
            'SERC_EGR': getattr(row, 'SERC_EGR', None),
            'fechaEgr': getattr(row, 'fechaEgr', None),
            'indexEgr': (getattr(row, 'fechaEgr', None)-min_fecha).days,
            'traslado': 0
        }

        first_tras = getattr(row, 'DIA_1_TRAS', None)
        if pd.notna(first_tras):
            new_row.update({
                'DIA_EGR': first_tras,
                'MES_EGR': getattr(row, 'MES_1_TRAS', None),
                'AÑO_EGR': getattr(row, 'ANO_1_TRAS', None),
                'AREAF_EGR': getattr(row, 'AREAF_1_TRAS', None),
                'SERC_EGR': getattr(row, 'SERC_1_TRAS', None),
                'fechaEgr': dt.datetime.strptime(
                    f"{int(getattr(row, 'ANO_1_TRAS'))}-{int(getattr(row, 'MES_1_TRAS'))}-{int(getattr(row, 'DIA_1_TRAS'))}",
                    "%Y-%m-%d"
                ) if pd.notna(getattr(row, 'ANO_1_TRAS', None)) else None,
                'traslado': 0
            })
            new_row.update({'indexEgr': (new_row.get('fechaEgr', None)-min_fecha).days})
            new_row.update({'estadia': new_row.get('indexEgr', None)-new_row.get('indexIng', None)})#+0.5
            new_row.update({'indice': l})
            rows.append(new_row.copy()) 
            k = 1
            for i, X in enumerate(steps):
                if pd.isna(getattr(row, f'DIA_{X}_TRAS', None)):
                    continue

                next_X = steps[i + 1] if i + 1 < len(steps) else None
                new_row = {col: getattr(row, col) for col in base_cols}

                new_row.update({
                    'DIA_ING': getattr(row, f'DIA_{X}_TRAS', None),
                    'MES_ING': getattr(row, f'MES_{X}_TRAS', None),
                    'AÑO_ING': getattr(row, f'ANO_{X}_TRAS', None),
                    'AREA_FUNC_I': getattr(row, f'AREAF_{X}_TRAS', None),
                    'UPC': int(getattr(row, f'AREAF_{X}_TRAS', None) in areasUPC_set),
                    'MB': int(getattr(row, f'AREAF_{X}_TRAS', None) in areasMB_set),
                    'SER_CLIN_I': getattr(row, f'SERC_{X}_TRAS', None),
                    'fechaIng': dt.datetime.strptime(
                        f"{int(getattr(row, f'ANO_{X}_TRAS'))}-{int(getattr(row, f'MES_{X}_TRAS'))}-{int(getattr(row, f'DIA_{X}_TRAS'))}",
                        "%Y-%m-%d"
                    ) if pd.notna(getattr(row, f'ANO_{X}_TRAS', None)) else None
                })

                new_row.update({'indexIng': (new_row.get('fechaIng', None)-min_fecha).days})
                new_row.update({'indice': l})

                if next_X and pd.notna(getattr(row, f'DIA_{next_X}_TRAS', None)):
                    new_row.update({
                        'DIA_EGR': getattr(row, f'DIA_{next_X}_TRAS', None),
                        'MES_EGR': getattr(row, f'MES_{next_X}_TRAS', None),
                        'AÑO_EGR': getattr(row, f'ANO_{next_X}_TRAS', None),
                        'AREAF_EGR': getattr(row, f'AREAF_{next_X}_TRAS', None),
                        'SERC_EGR': getattr(row, f'SERC_{next_X}_TRAS', None),
                        'fechaEgr': dt.datetime.strptime(
                            f"{int(getattr(row, f'ANO_{next_X}_TRAS'))}-{int(getattr(row, f'MES_{next_X}_TRAS'))}-{int(getattr(row, f'DIA_{next_X}_TRAS'))}",
                            "%Y-%m-%d"
                        ) if pd.notna(getattr(row, f'ANO_{next_X}_TRAS', None)) else None,
                        'traslado': k
                    })
                    k +=1
                    new_row.update({'indexEgr': (new_row.get('fechaEgr', None)-min_fecha).days})
                    new_row.update({'estadia': new_row.get('indexEgr', None)-new_row.get('indexIng', None)})

                else:
                    new_row.update(last_egr)
                    new_row.update({'estadia': new_row.get('indexEgr', None)-new_row.get('indexIng', None)}) #+0.5
                    new_row.update({
                        'traslado': k
                    })

                rows.append(new_row.copy())

        else:
            new_row.update(last_egr)
            new_row.update({'estadia': new_row.get('indexEgr', None)-new_row.get('indexIng', None)})
            rows.append(new_row.copy())

        l+=1


    return pd.DataFrame(rows)


###############################################################################
###############################################################################
###############################################################################
############################### PRE TRATAMIENTO ###############################
###############################################################################
###############################################################################
###############################################################################


###############################################################################
############################### Atenciones de Urgencias #######################
###############################################################################


# dics_urgencias = []
# for d in cf.dic_urgencias:
#     dics_urgencias.append((pd.read_excel(folder_path_input+d, sheet_name='Anexo 1', skiprows=5))[['ID CAUSA', 'GLOSACAUSA']])
# dics_urgencias = pd.concat(dics_urgencias).drop_duplicates()
# dics_urgencias.to_csv(f'{folder_path}/dics_urgencias.csv')
# urgencias = []
# col_ur = ['IdEstablecimiento','NEstablecimiento','IdCausa','GlosaCausa',
#                'Total','Menores_1','De_1_a_4','De_5_a_14','De_15_a_64','De_65_y_mas',
#                 'fecha','semana','GLOSATIPOESTABLECIMIENTO','GLOSATIPOATENCION','GlosaTipoCampana']
# for u in cf.urgencias:
#     df_urgencias = (pd.read_csv(folder_path_input+u,encoding= 'latin1', sep = ';', low_memory=False))[col_ur]
#     df_urgencias = df_urgencias[df_urgencias['IdCausa'].isin(cf.urgencia_causas)]
#     urgencias.append(df_urgencias)

# urgencias = pd.concat(urgencias)
# urgencias = urgencias.assign(categoria_ICD10 = lambda urgencias: urgencias['IdCausa'].map(cf.urgencias_causas_dic))
# urgencias.to_parquet(f'{folder_path}/REM_urgencias/urgencias_respiratory_causes.parquet', engine='pyarrow', compression='snappy')


# # Final y tabla resumen
# urgencias['fecha'] = pd.to_datetime(urgencias['fecha'], format='%d/%m/%Y')
# urgencias = (urgencias
#              .assign(year = lambda urgencias: urgencias['fecha'].dt.year,
#                      month = lambda urgencias: urgencias['fecha'].dt.month,
#                      week = lambda urgencias: urgencias['fecha'].dt.isocalendar().week))
# urgencias_c = (urgencias
#                   .sort_values(['year','month'])
#                   .groupby(['categoria_ICD10','year','month'])
#                   .agg({'Total':'sum',
#                         'Menores_1':'sum',
#                         'De_1_a_4':'sum'})   #,'De_65_y_mas':'sum'})
#                   .reset_index()
#                   .drop_duplicates())
# outpatient_table = (urgencias_c
#                     .sort_values(['year'])
#                     .groupby(['categoria_ICD10','year'])
#                     .agg({'Menores_1':'sum',
#                           'De_1_a_4':'sum'})   
#                     .reset_index()
#                     .drop_duplicates()
#                     .assign(Total = lambda urgencias_c: urgencias_c['Menores_1']+urgencias_c['De_1_a_4']))
# outpatient_table_icd10 = pd.pivot_table(outpatient_table, values='Total', index=['categoria_ICD10'], columns=['year'], aggfunc=np.sum, fill_value=0)
# outpatient_table_age = (outpatient_table[['year','Menores_1','De_1_a_4']]
#                         .groupby(['year'])
#                         .agg({'Menores_1':'sum',
#                               'De_1_a_4':'sum'})   
#                         .reset_index()
#                         .drop_duplicates()
#                         .set_index('year')).T
# outpatient_table_icd10.to_excel(folder_path+'/Paper/outpatient_table_icd10.xlsx')
# outpatient_table_age.to_excel(folder_path+'/Paper/outpatient_table_age.xlsx')

###############################################################################
################################# Nacimientos #################################
###############################################################################

path_actual = Path.cwd()
path_data = path_actual.parent/'Data' 

nac_24 = pd.read_csv(path_data / 'COHORTE_NIRSE_ACTUALIZADA_23_06_2025_ENCR.csv',encoding= 'latin1', sep = ';', low_memory=False)
nac = pd.read_csv(path_data / 'NAC_DEF_EGRE3.csv',encoding= 'latin1', sep = ';', low_memory=False)
nac_24 = (nac_24
            .assign(fechaIng=lambda x: pd.to_datetime(x['FECHA_INGRESO'],dayfirst=True),
                    fechaNac=lambda x: pd.to_datetime(x['FECHA_NACIMIENTO'],dayfirst=True),
                    fechaEgr=lambda x: pd.to_datetime(x['FECHA_EGRESO'],dayfirst=True),
                    fechaDef=lambda x: pd.to_datetime(x['FECHA_DEFUNCION'],dayfirst=True),
                    fechaInm=lambda x: pd.to_datetime(x['FECHA_INMUNIZACION'],dayfirst=True))
            
            .assign(DIA_ING=lambda x: x['fechaIng'].dt.day,
                    MES_ING=lambda x: x['fechaIng'].dt.month,
                    ANO_ING=lambda x: x['fechaIng'].dt.year,
                    DIA_EGR=lambda x: x['fechaEgr'].dt.day,
                    MES_EGR=lambda x: x['fechaEgr'].dt.month,
                    ANO_EGR=lambda x: x['fechaEgr'].dt.year,
                    D_NAC=lambda x: x['fechaNac'].dt.day,
                    M_NAC=lambda x: x['fechaNac'].dt.month,
                    A_NAC=lambda x: x['fechaNac'].dt.year))
nac = (nac
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

nac.to_parquet(f'{folder_path}/Nacimientos/nacimientos.parquet', engine='pyarrow', compression='snappy')
nac_24.to_parquet(f'{folder_path}/Nacimientos/nacimientos_2024.parquet', engine='pyarrow', compression='snappy')
nacimientos = pd.concat([nac,nac_24],ignore_index=True)
# Nacimientos con las fechas del estudio
nacimientos = nacimientos[ 
    (( nacimientos['A_NAC']<2023 ) 
    |
    ( (nacimientos['A_NAC']==2023) & (nacimientos['M_NAC']<11) )
    |
    ( (nacimientos['A_NAC']==2023) & (nacimientos['M_NAC']==11) & (nacimientos['D_NAC']<24) )  
    )]
#nacimientos.to_parquet(f'{folder_path}/Paper/IMPACTO/nacimientos.parquet', engine='pyarrow', compression='snappy')

# En nacimientos no se llego a los mismos valores reportados por el paper
# (no estaban todos los filtros que se habían ocupado en el código entregado)
# por lo que filtre para que se pareciera al menos a los del INE

# Filtro 1
nacimientos = (nacimientos.dropna(subset=['SEMANAS','PESO']))
nacimientos= (nacimientos.assign(year_nac = lambda x: x['fechaNac'].dt.year))

# Filtro 2: Se estudio la distribución de peso y semanas para filtrar
dist_semanas_o = (nacimientos
                .groupby(['SEMANAS','year_nac'])
                .agg({'RUN':'count'})
                .reset_index()
                .assign(posible = lambda x: (x['SEMANAS']>=24)&(x['SEMANAS']<=42)))
dist_peso_o = (nacimientos
                .groupby(['PESO','year_nac'])
                .agg({'RUN':'count'})
                .reset_index()
                .assign(posible = lambda x: (x['PESO']>=1500)&(x['PESO']<=8000)))
nacimientos2 = (nacimientos
                [(nacimientos['PESO']>0)&(nacimientos['PESO']<6000)&
                 (nacimientos['SEMANAS']>4)&(nacimientos['SEMANAS']<60)] )

# Filtro 3: Se filtro a los nacidos que fallecen el mismo día y los ruts que estaban duplicados
run_fallecidos = nacimientos2[(nacimientos2['fechaDef']==nacimientos2['fechaNac'])].RUN.to_list()
nacimientos3 = nacimientos2.query('~RUN.isin(@run_fallecidos)')
nacimientos_run = nacimientos3.groupby('RUN').size().reset_index(name='count').sort_values('count',ascending=False)
run_ok = nacimientos_run.query('count==1').RUN.to_list()
run_not_ok = nacimientos_run.query('count>1').RUN.to_list()
nacimientos_ok = nacimientos3.query('RUN.isin(@run_ok)')
nacimientos_nok = nacimientos3.query('RUN.isin(@run_not_ok)')
nacimientos_ok2 = (nacimientos_nok
                   .drop_duplicates(subset=['RUN','fechaNac','SEMANAS','PESO','TALLA'],keep='first'))
temp = nacimientos_ok2.groupby('RUN').size().reset_index(name='count').sort_values('count',ascending=False)
temp_ok = temp.query('count==1').RUN.to_list()
temp_nok = temp.query('count>1').RUN.to_list()
nacimientos_nok2 = nacimientos_ok2[nacimientos_ok2['RUN'].isin(temp_nok)]
nacimientos_ok2 = nacimientos_ok2[nacimientos_ok2['RUN'].isin(temp_ok)]

# Final y tabla resumen
nacimientos_final = pd.concat([nacimientos_ok,nacimientos_ok2])
nactual = (nacimientos_final
           .groupby('year_nac')
           .agg({'RUN':'count'})
           .reset_index()
           .rename(columns={'RUN':'total_actual'})
           .query('year_nac.isin([2019,2020,2021,2022,2023])')
           .sort_values('year_nac',ascending=True))
npaper = [245440, 217983, 197219, 208095, 173920]
nine = [210188, 194978, 177273, 189310, 171992]
nactual['total_paper'] = npaper
nactual['total_ine'] = nine
nactual = (nactual
           .assign(error = lambda x: round(((x['total_ine']-x['total_actual'])/x['total_ine'])*100,2)))
nacimientos_res = nactual.set_index('year_nac').T
nacimientos_res.to_excel(path_data/'nacimientos_contabilizados_amal.xlsx')

###############################################################################
################################### EGRESOS ###################################
###############################################################################
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
        .assign(in_season = lambda df: inseason(df['M_NAC'],df['A_NAC'],df['ANO_ING'])))
comunas_dic = pd.read_csv(f'{folder_path}/comunas_dic.csv').rename(columns={'comuna_str':'comuna_paciente_str','Region':'Region_paciente'})
df = (df.merge(comunas_dic,how='left',on='COMUNA'))
df = df.drop(columns=['Unnamed: 0'])
# egresos_full.parquet
df.to_parquet(folder_path+'/egresos_all.parquet', engine='pyarrow', compression='snappy')
#egresos_5.parquet
df_5 = df[df['dias_vida']<=365*5+3]
df_5.to_parquet(folder_path+'/egresos_5.parquet', engine='pyarrow', compression='snappy')

min_fecha = cf.min_fecha
max_fecha = cf.max_fecha
season = []
years = []
years_ = [2019,2020,2021,2022,2023]
for year in [2019,2020,2021,2022,2023]:
    dates1 = pd.to_datetime({'year': [year], 'month': [4], 'day': [1]}, format='%Y-%m-%d')[0] # fecha inicio temporada 01/04
    dates2 = pd.to_datetime({'year': [year], 'month': [9], 'day': [30]}, format='%Y-%m-%d')[0] # fecha fin temporada 30/09
    datey1 = pd.to_datetime({'year': [year], 'month': [1], 'day': [1]}, format='%Y-%m-%d')[0] 
    datey2 = pd.to_datetime({'year': [year], 'month': [12], 'day': [31]}, format='%Y-%m-%d')[0] 
    season.append(((dates1 - min_fecha).days, (dates2 - min_fecha).days))
    years.append(((datey1 - min_fecha).days, (datey2 - min_fecha).days))
# Tabla del paper
years_nirse = [2019,2020,2021,2022,2023]
menor_1 = [79357,47790,46451,60977,45398]
de_1_a_4 = [11222,12820,23250,39419,28088]
male = [51032,33689,39294,57224,40886]
female = [39325,26902,30388,43141,32581]
other = [16,19,19,31,19]
influenza = [914,35,56,808,538]
pneumonia = [10763,825,2137,10596,13065]
upper_respiratory_tract_infection = [2783,1099,1482,3616,2292]
obstructive_pulmonary_disease = [1802,392,892,2909,2982]
acute_bronchitis_bronchiolitis = [14707,1136,5067,13325,13437]
covid_19 = [0,800,693,2487,684]
other_respiratory_conditions = [4337,1168,3098,11006,12316]
non_respiratory_disease = [60723,55780,58495,64970,40504]
df_paper = pd.DataFrame({'year':years_nirse,
                            'menor_1':menor_1,
                            'de_1_a_4':de_1_a_4,
                            'Male':male,
                            'Female':female,
                            'desconocido':other,
                            'influenza':influenza,
                            'pneumonia':pneumonia,
                            'upper_respiratory_tract_infection':upper_respiratory_tract_infection,
                            'obstructive_pulmonary_disease':obstructive_pulmonary_disease,
                            'acute_bronchitis_bronchiolitis':acute_bronchitis_bronchiolitis,
                            'covid':covid_19,
                            'other_respiratory_conditions':other_respiratory_conditions,
                            'non_respiratory_disease':non_respiratory_disease})
df_paper.set_index('year').T.to_excel(folder_path+'/Paper/EGRESOS/tabla1_paper.xlsx')
# Tabla egresos generada por réplica código gonzalo y por amal
df = (pd.read_parquet(f'{folder_path}/egresos_5.parquet'))
nacimientos = pd.read_parquet(f'{folder_path}/Nacimientos/nacimientos.parquet')
cols_nac = ['RUN','fechaNac','fechaDef','SEMANAS','PESO',
            'TALLA','D_NAC','M_NAC','A_NAC']
nacimientos2 = nacimientos[cols_nac]
nacimientos2 = nacimientos2.drop_duplicates(subset=['RUN'])
df = df.merge(nacimientos2,on=['RUN','fechaNac','A_NAC','D_NAC','M_NAC'],how='left')
# Temporadas
fechas_inicio = {y: pd.Timestamp(y, 4, 1) for y in range(2018, 2024)}
fechas_fin = {y: pd.Timestamp(y, 9, 30) for y in range(2019, 2024)}
df = (
    df
    .assign(riesgo=lambda x: (x['SEMANAS'] < 32) | (x['PESO'] < 1500))
)
for y in range(2019, 2024):
    col_catchup = f'catchup{y}'
    col_posible = f'posible_{y}'
    col_riesgo = f'riesgo_{y}'
    col_inseason = f'inseason_{y}'

    df[col_posible] = (
        (df['fechaNac'] >= fechas_inicio[y - 1]) &
        (df['fechaNac'] <= fechas_fin[y])
    )
    df[col_inseason] = (
        (df['fechaNac'] >= fechas_inicio[y]) &
        (df['fechaNac'] <= fechas_fin[y])
    )
    df[col_riesgo] = (
        df[col_posible] &
        df['riesgo']
    )
    df[col_catchup] = (
        (df['fechaNac'] >= fechas_inicio[y - 1]) &
        (df['fechaNac'] < fechas_inicio[y])
    )
# Filtros
df = df[(df['year'] > 2018) & 
        (df['year'] <=2023) &
        (~df['ESTAB'].isna()) & 
        (df['dias_vida'] >= 0) & 
        (df['DIAS_ESTAD'] >= 0)]

# Calculo de edad según código gonzalo
mask = df['fechaIng'].notna() & df['fechaNac'].notna()
df['edad_meses'] = np.nan
df_valid = df.loc[mask].copy()
rel_diff = df_valid.apply(lambda row: relativedelta(row['fechaIng'], row['fechaNac']), axis=1)
df.loc[mask, 'edad_meses'] = rel_diff.map(lambda rd: rd.years * 12 + rd.months)
df = df[df['edad_meses']<=48]

# Diagnosticos según código gonzalo
diag_col = [col for col in df.columns if col.startswith("DIAG")]
rsv_0_set = set(cf.rsv_0)
rsv_1_set = set(cf.rsv_1)
rsv_2_set = set(cf.rsv_2)
df['vrs0'] = df[diag_col].isin(rsv_0_set).any(axis=1)
df['vrs1'] = df[diag_col].isin(rsv_1_set).any(axis=1)
df['vrs2'] = df[diag_col].isin(rsv_2_set).any(axis=1)
diag_codes = pd.read_excel(f'{folder_path_input}/IMPACTO/codigos.xlsx', sheet_name='Diag2')
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
irag = set(pd.read_csv(f'{folder_path_input}/IMPACTO/codigos_IRAG.csv').IRAG.to_list())
df['irag'] = df[diag_col].isin(irag).any(axis=1)

# Clasificación de edad en años
df['EDAD_ANOS'] = np.select(
    [
        df['edad_meses'] < 12,
        (df['edad_meses'] >= 12) & (df['edad_meses'] < 24),
        (df['edad_meses'] >= 24) & (df['edad_meses'] < 36),
        (df['edad_meses'] >= 36) & (df['edad_meses'] < 48),
        (df['edad_meses'] >= 48) & (df['edad_meses'] < 60),
        df['edad_meses'] >= 60
    ],
    [0, 1, 2, 3, 4, 5]
)
df = (df.assign(menores_1 = lambda x: 1*(x['edad_meses']<=12),
                de_1_a_4 = lambda x: 1*((x['edad_meses']>12)&(x['edad_meses']<=48))))
# Filtros de fechas
fecha_ing_min = pd.to_datetime({'year': [2019], 'month': [1], 'day': [1]}, format='%Y-%m-%d')[0]
df = df[df['fechaIng']>fecha_ing_min]
fecha_ing_max = pd.to_datetime({'year': [2023], 'month': [9], 'day': [30]}, format='%Y-%m-%d')[0]
df = df[df['fechaIng']<=fecha_ing_max]

################## Datos Egresos Formato Gonzalo ##################
diseases = df.copy()
# Forma de filtrar duplicados
diseases['ing'] = diseases.groupby(['RUN', 'fechaIng']).ngroup()
errors = diseases['ing'].value_counts()
errors = errors[errors > 1].index.tolist()
diseases = diseases[~diseases['ing'].isin(errors)]
diseases = diseases[diseases['fechaIng'].notna()].copy()
areaf_cols = diseases.filter(like='AREAF').columns
diseases['UPC_MB'] = diseases[areaf_cols].isin(set(cf.areasUPC+cf.areasMB)).any(axis=1)
fecha_nac_min = pd.to_datetime({'year': [2018], 'month': [1], 'day': [1]}, format='%Y-%m-%d')[0]
diseases = diseases[diseases['fechaNac']>fecha_nac_min]
diseases = diseases[diseases['UPC_MB'] & diseases['fechaIng'].notna()].copy()
diseases['ageGroup'] = diseases['edad_meses'].apply(lambda x: 'menor_1' if x <= 12 else ('de_1_a_4' if x <= 48 else '4+'))
# Conteo por año y grupo de edad
age_group_summary = diseases.groupby(['ANO_ING', 'ageGroup']).size().reset_index(name='n')
age_group_summary_pivot = age_group_summary.pivot(index='ageGroup', columns='ANO_ING', values='n')
# Conteo por sexo y año (menores de 48 meses)
sex_summary = (diseases[diseases['UPC_MB'] & (diseases['edad_meses'] <= 48)]
               .replace(['intersex','desconocido'],'desconocido')
               .groupby(['ANO_ING', 'sexo_str'])
               .size()
               .reset_index(name='n'))
sex_summary_pivot = sex_summary.pivot(index='sexo_str', columns='ANO_ING', values='n')
# Conteo de enfermedades por año (solo menores de 48 meses en UPC)
disease_year_summary = (diseases[diseases['UPC_MB'] & (diseases['edad_meses'] <= 48)]
                        .groupby('ANO_ING')[disease_cols + ['non_respiratory_disease']].sum().T)
formato_g = pd.concat([pd.concat([age_group_summary_pivot,sex_summary_pivot,disease_year_summary])])
estab_publicos = [
    set(diseases[(diseases['origen_estab'] == 'public') & (diseases['ANO_ING'] == year)].ESTAB.unique()) 
    for year in [2019, 2020, 2021, 2022, 2023]
]
estab_pub_u = set.union(*estab_publicos)
estab_pub = set.intersection(*estab_publicos)
estab_privados = [
    set(diseases[(diseases['origen_estab'] == 'private') & (diseases['ANO_ING'] == year)].ESTAB.unique()) 
    for year in [2019, 2020, 2021, 2022, 2023]
]
estab_priv_q = [len(e) for e in estab_privados]
estab_priv_u = set.union(*estab_privados)
print(f'Establecimientos públicos: {len(estab_pub_u)} vs 156 \nEstablecimientos privados: {len(estab_priv_u)} vs 78')
formato_g.to_excel(folder_path+'/Paper/EGRESOS/tabla1_gonzalo.xlsx')
diseases.to_parquet(folder_path+'/Paper/IMPACTO/egresos_gonzalo.parquet')


################## Datos Egresos Formato Amal ##################
diseases_amal = df.copy()
diseases_amal = diseases_amal.drop_duplicates(['RUN','fechaIng'],keep='first')
areaf_cols = diseases_amal.filter(like='AREAF').columns
diseases_amal['UPC_MB'] = diseases_amal[areaf_cols].isin(set(cf.areasUPC+cf.areasMB)).any(axis=1)
diseases_amal = diseases_amal[diseases_amal['UPC_MB'] & diseases_amal['fechaIng'].notna()].copy()
diseases_amal['ageGroup'] = diseases_amal['edad_meses'].apply(lambda x: 'menor_1' if x <= 12 else ('de_1_a_4' if x <= 48 else '4+'))
# Conteo por año y grupo de edad
age_group_summary = diseases_amal.groupby(['ANO_ING', 'ageGroup']).size().reset_index(name='n')
age_group_summary_pivot_amal = age_group_summary.pivot(index='ageGroup', columns='ANO_ING', values='n')
# Conteo por sexo y año (menores de 48 meses)
sex_summary = (diseases_amal[diseases_amal['UPC_MB'] & (diseases_amal['edad_meses'] <= 48)]
               .replace(['intersex','desconocido'],'desconocido')
               .groupby(['ANO_ING', 'sexo_str'])
               .size()
               .reset_index(name='n'))
sex_summary_pivot_amal = sex_summary.pivot(index='sexo_str', columns='ANO_ING', values='n')
# Conteo de enfermedades por año (solo menores de 48 meses en UPC)
disease_year_summary_amal = (diseases_amal[diseases_amal['UPC_MB'] & (diseases_amal['edad_meses'] <= 48)]
                        .groupby('ANO_ING')[disease_cols + ['non_respiratory_disease']].sum().T)
formato_a = pd.concat([pd.concat([age_group_summary_pivot_amal,sex_summary_pivot_amal,disease_year_summary_amal])])
estab_publicos = [
    set(diseases_amal[(diseases_amal['origen_estab'] == 'public') & (diseases_amal['ANO_ING'] == year)].ESTAB.unique()) 
    for year in [2019, 2020, 2021, 2022, 2023]
]
estab_pub_u = set.union(*estab_publicos)
estab_pub = set.intersection(*estab_publicos)
estab_privados = [
    set(diseases_amal[(diseases_amal['origen_estab'] == 'private') & (diseases_amal['ANO_ING'] == year)].ESTAB.unique()) 
    for year in [2019, 2020, 2021, 2022, 2023]
]
estab_priv_q = [len(e) for e in estab_privados]
estab_priv_u = set.union(*estab_privados)
print(f'Establecimientos públicos: {len(estab_pub_u)} vs 156 \nEstablecimientos privados: {len(estab_priv_u)} vs 78')
formato_a.to_excel(folder_path+'/Paper/EGRESOS/tabla1_amal_de14.xlsx')
diseases_amal.to_parquet(folder_path+'/Paper/IMPACTO/egresos_amal_04.parquet')
## Datos Formato Amal elegibles por año
diseases_amal2 = df.copy()
diseases_amal2 = diseases_amal2.drop_duplicates(['RUN','fechaIng'],keep='first')
areaf_cols = diseases_amal2.filter(like='AREAF').columns
diseases_amal2['UPC_MB'] = diseases_amal2[areaf_cols].isin(set(cf.areasUPC+cf.areasMB)).any(axis=1)
diseases_amal2 = diseases_amal2[diseases_amal2['UPC_MB'] & diseases_amal2['fechaIng'].notna()].copy()
diseases_amal2['ageGroup'] = diseases_amal2['edad_meses'].apply(lambda x: 'menor_1' if x <= 12 else ('de_1_a_2' if x <= 24 else '2+'))
col_posible = [f'posible_{y}' for y in range(2019, 2024)]
# Conteo por año y grupo de edad
age_group_summary2 = diseases_amal2.groupby(['ANO_ING', 'ageGroup']+col_posible).size().reset_index(name='n')
# Conteo por sexo y año (menores de 48 meses)
sex_summary2 = (diseases_amal2[diseases_amal2['UPC_MB'] & (diseases_amal2['edad_meses'] <= 48)]
               .replace(['intersex','desconocido'],'desconocido')
               .groupby(['ANO_ING', 'sexo_str']+col_posible)
               .size()
               .reset_index(name='n'))
# Conteo de enfermedades por año (solo menores de 48 meses en UPC)
disease_year_summary_amal2 = (diseases_amal2[diseases_amal2['UPC_MB'] & (diseases_amal2['edad_meses'] <= 48)]
                        .groupby(['ANO_ING']+col_posible)[disease_cols + ['non_respiratory_disease']].sum().reset_index())
age_group_summary2.to_excel(folder_path+'/Paper/EGRESOS/tabla1_amal_elegibles_ageGroup.xlsx')
sex_summary2.to_excel(folder_path+'/Paper/EGRESOS/tabla1_amal_elegibles_sex.xlsx')
disease_year_summary_amal2.to_excel(folder_path+'/Paper/EGRESOS/tabla1_amal_elegibles_diseases.xlsx')
estab_publicos = [
    set(diseases_amal2[(diseases_amal2['origen_estab'] == 'public') & (diseases_amal2['ANO_ING'] == year)].ESTAB.unique()) 
    for year in [2019, 2020, 2021, 2022, 2023]
]
estab_pub_u = set.union(*estab_publicos)
estab_pub = set.intersection(*estab_publicos)
estab_privados = [
    set(diseases_amal2[(diseases_amal2['origen_estab'] == 'private') & (diseases_amal2['ANO_ING'] == year)].ESTAB.unique()) 
    for year in [2019, 2020, 2021, 2022, 2023]
]
estab_priv_q = [len(e) for e in estab_privados]
estab_priv_u = set.union(*estab_privados)
print(f'Establecimientos públicos: {len(estab_pub_u)} vs 156 \nEstablecimientos privados: {len(estab_priv_u)} vs 78')
diseases_amal2.to_parquet(folder_path+'/Paper/IMPACTO/egresos_amal_posibles.parquet')

###############################################################################
###############################################################################
###############################################################################
################################### CURVAS  ###################################
###############################################################################
###############################################################################
###############################################################################

###############################################################################
################################# Nacimientos #################################
###############################################################################
# Fechas de inicio y fin de temporada 
date_inicio = {year: pd.Timestamp(year, 4, 1) for year in range(2018, 2024)}
date_fin = {year: pd.Timestamp(year, 9, 30) for year in range(2019, 2024)}
# Columnas de temporada y riesgo
def agregar_temporadas(df):
    df = df.copy()
    df['riesgo'] = (df['SEMANAS'] < 32) | (df['PESO'] < 1500)
    for year in range(2019, 2024):
        df[f'posible_{year}'] = (df['fechaNac'] >= date_inicio[year - 1]) & (df['fechaNac'] <= date_fin[year])
        df[f'inseason_{year}'] = (df['fechaNac'] >= date_inicio[year]) & (df['fechaNac'] <= date_fin[year])
        df[f'riesgo_{year}'] = df[f'posible_{year}'] & df['riesgo']
    return df
nacimientos_final_riesgo = agregar_temporadas(nacimientos_final)
print(f'Pacientes riesgo temporada 2023: {nacimientos_final_riesgo.riesgo_2023.sum()}')
print(f'Posibles riesgosos temporada 2023 (menos de un año al comienzo): {nacimientos_final_riesgo.posible_2023.sum()}')
print(f'Porcentaje 2023: {round(nacimientos_final_riesgo.riesgo_2023.sum()/nacimientos_final_riesgo.posible_2023.sum()*100,2)}%')
print(f'Nacimientos en la temporada 2023 (inseason): ' ,nacimientos_final_riesgo.query('riesgo_2023==0').inseason_2023.sum())

#### Nacimientos catchup estrategía
def catchup(fechaNac, y):
    fechas = [pd.Timestamp(year=y, month=m, day=1) for m in range(4, 0, -1)] + \
             [pd.Timestamp(year=y-1, month=m, day=1) for m in range(12, 3, -1)]
    for i, fecha_limite in enumerate(fechas):
        if fechaNac >= fecha_limite:
            return i
    return -1
for año in [2019, 2022, 2023]:
    col_catchup = f'catchup{año}'
    col_posible = f'posible_{año}'
    col_inseason = f'inseason_{año}'
    nacimientos_final_riesgo[col_catchup] = nacimientos_final_riesgo.apply(
        lambda row: catchup(row['fechaNac'], año) if row[col_posible] and not row[col_inseason] else -1,
        axis=1
    )
nacimientos_final_riesgo.to_parquet(f'{folder_path}/Paper/FINALES/Tablas/nacimientos.parquet')

### Curva de Nacimientos
curva_nacimientos = pd.DataFrame({'year':[2019,2022,2023]})
curva_nacimientos['riesgo_inseason'] = [
        nacimientos_final_riesgo.query(f'{col_inseason}==1')[col_riesgo].sum()
        for (col_riesgo,col_inseason)
        in [('riesgo_2019','inseason_2019'),('riesgo_2022','inseason_2022'),('riesgo_2023','inseason_2023')]]

curva_nacimientos['riesgo_catchup'] = [
        nacimientos_final_riesgo.query(f'{col_inseason}==0')[col_riesgo].sum()
        for (col_riesgo,col_inseason)
        in [('riesgo_2019','inseason_2019'),('riesgo_2022','inseason_2022'),('riesgo_2023','inseason_2023')]]

curva_nacimientos['inseason'] = [
        nacimientos_final_riesgo.query(f'{col_riesgo}==0')[col_inseason].sum()
        for (col_riesgo,col_inseason)
        in [('riesgo_2019','inseason_2019'),('riesgo_2022','inseason_2022'),('riesgo_2023','inseason_2023')] ]

curva_nacimientos['posibles'] = [
        nacimientos_final_riesgo.query(f'{col_posible}==0')[col_inseason].sum()
        for col_posible
        in ['posible_2019','posible_2022','posible_2023'] ]

for k in range(0,13):
    col = f'catchup{k}'
    curva_nacimientos[col] = [
            (nacimientos_final_riesgo
            .query(f'{col_riesgo}==0')
            .query(f'{col_posible}==1')
            .query(f'{col_inseason}==0')
            .query(f'{col_catchup}>=0')
            .query(f'{col_catchup}<={k}'))[col_catchup].size
            for (col_riesgo,col_inseason,col_catchup,col_posible)
            in [('riesgo_2019','inseason_2019','catchup2019','posible_2019'),
                ('riesgo_2022','inseason_2022','catchup2022','posible_2022'),
                ('riesgo_2023','inseason_2023','catchup2023','posible_2023')]]
curva_nacimientos.to_excel(f'{folder_path}/Paper/FINALES/Tablas/curva_nacimientos.xlsx')

###############################################################################
################################### EGRESOS ###################################
###############################################################################

# Solo se miran egresos LRTI
# Se crean curvas de MB y UPC con el formato gonzalo y amal
def curvaHosp_v2(dataframe):
    
    max_ind = dataframe['indexEgr'].max()
    diag = 1


    df_ing_UPC = (dataframe[dataframe['UPC']==1]
                            .groupby(['indexIng','vrs1'])
                            .agg({'RUN':'count'})
                            .reset_index()
                            .rename(columns={'indexIng':'index', 'RUN': 'Count'})
                            .drop_duplicates())
    df_egr_UPC = (dataframe[dataframe['UPC']==1]
                            .groupby(['indexEgr','vrs1'])
                            .agg({'RUN':'count'})
                            .reset_index()
                            .rename(columns={'indexEgr':'index', 'RUN': 'Count'})
                            .drop_duplicates())
    df_ing_MB = (dataframe[dataframe['MB']==1]
                            .groupby(['indexIng','vrs1'])
                            .agg({'RUN':'count'})
                            .reset_index()
                            .rename(columns={'indexIng':'index', 'RUN': 'Count'})
                            .drop_duplicates())
    df_egr_MB = (dataframe[dataframe['MB']==1]
                            .groupby(['indexEgr','vrs1'])
                            .agg({'RUN':'count'})
                            .reset_index()
                            .rename(columns={'indexEgr':'index', 'RUN': 'Count'})
                            .drop_duplicates())
    
    df_egr_UPC['Count'] = -df_egr_UPC['Count']
    df_egr_MB['Count'] = -df_egr_MB['Count']
    df_mov_UPC = pd.concat([df_ing_UPC,df_egr_UPC]).sort_values('index',ascending=True)
    df_mov_UPC['UPC'] = 1
    df_mov_UPC['MB'] = 0
    df_mov_MB = pd.concat([df_ing_MB,df_egr_MB]).sort_values('index',ascending=True)
    df_mov_MB['MB'] = 1
    df_mov_MB['UPC'] = 0
    df_mov = pd.concat([df_mov_MB,df_mov_UPC]).sort_values('index',ascending=True)
    

    curvas = np.zeros((max_ind+1, diag+1, 2)) 
    # fecha (min_ind,_max_ind),diagnostico (diccionario mapper), upc (1 critico / 0 media basica)

    k = np.zeros((diag+1,2))
    i = 0

    for row in df_mov.itertuples(index=False):
        d = getattr(row,'vrs1',None)*1
        upc = getattr(row,'UPC',None)

        if i < getattr(row,'index',None):
            curvas[i:getattr(row,'index',None)+1,:,:] = k
        
        i = getattr(row,'index',None)
        curvas[i,d,upc] += getattr(row,'Count',None)
        k = curvas[i]
    
    return curvas

############################ Egresos LRTI Gonzalo ############################

egresos_lrti = diseases[(diseases['irag'])|(diseases['vrs1'])]
for año in [2019, 2022, 2023]:
    col_catchup = f'catchup2{año}'
    col_posible = f'posible_{año}'
    col_inseason = f'inseason_{año}'
    egresos_lrti[col_catchup] = egresos_lrti.apply(
        lambda row: catchup(row['fechaNac'], año) if row[col_posible] and not row[col_inseason] else -1,
        axis=1
    )
### Egresos elegibles estrategía
for y in [2019,2022,2023]:
    for k in range(0,13):
        col = f'elegible{k}_{y}'
        col_catchup = f'catchup2{y}'
        col_posible = f'posible_{y}'
        col_inseason = f'inseason_{y}'
        col_riesgo = f'riesgo_{y}'
        egresos_lrti[col] =  ((egresos_lrti[col_riesgo]) | (egresos_lrti[col_inseason]) |
                               ((egresos_lrti[col_catchup]>=0) & (egresos_lrti[col_catchup]<=k)))
                
egresos_lrti_estadias = transform_translations(egresos_lrti,cf.areasUPC,cf.areasMB,min_fecha)
egresos_lrti_estadias_2019 = egresos_lrti_estadias.query('indexIng>=@season[0][0]').query('indexIng<=@season[0][1]')
egresos_lrti_estadias_2022 = egresos_lrti_estadias.query('indexIng>=@season[3][0]').query('indexIng<=@season[3][1]')
egresos_lrti_estadias_2023 = egresos_lrti_estadias.query('indexIng>=@season[4][0]').query('indexIng<=@season[4][1]')

def construir_df(curvas, elegible, cama):
    df = pd.DataFrame([
        c for c in curvas[cama]['elegible' if elegible else 'no_elegible']
    ], columns=['estrategia', 'year', 'curva'])

    df_expandido = df['curva'].apply(pd.Series)
    df_final = pd.concat([df.drop(columns='curva'), df_expandido], axis=1)
    df_final['cama'] = cama
    df_final['elegible'] = int(elegible)
    df_final['VRS'] = 1
    return df_final
# Configuración general
años = [2019, 2022, 2023]
campos = {'UPC': 1, 'MB': 0}
rango_k = range(13)
season_map = {2019: 0, 2022: 3, 2023: 4}
egresos_map = {
    2019: egresos_lrti_estadias_2019,
    2022: egresos_lrti_estadias_2022,
    2023: egresos_lrti_estadias_2023,
}
### UPC VRS
curvas_upc = {cama: {'elegible': [], 'no_elegible': []} for cama in campos}
for cama, idx_cama in campos.items():
    for year in años:
        df = egresos_map[year]
        s_ini, s_fin = season[season_map[year]]

        for k in rango_k:
            col = f'elegible{k}_{year}'

            curva_elegible = curvaHosp_v2(df.query(f'{col}==True'))[s_ini:s_fin, 1, idx_cama]
            curva_no_elegible = curvaHosp_v2(df.query(f'{col}==False'))[s_ini:s_fin, 1, idx_cama]

            curvas_upc[cama]['elegible'].append((k, year, curva_elegible))
            curvas_upc[cama]['no_elegible'].append((k, year, curva_no_elegible))
df_lrti_vrs_UPC = pd.concat([
    construir_df(curvas_upc, elegible=True, cama='UPC'),
    construir_df(curvas_upc, elegible=False, cama='UPC')
])
df_lrti_vrs_UPC.to_csv(f'{folder_path}/Paper/FINALES/df_lrti_vrs_UPC_gonzalo.csv')
### MB VRS
curvas_mb = {cama: {'elegible': [], 'no_elegible': []} for cama in campos}
for cama, idx_cama in campos.items():
    for year in años:
        df = egresos_map[year]
        s_ini, s_fin = season[season_map[year]]

        for k in rango_k:
            col = f'elegible{k}_{year}'

            curva_elegible = curvaHosp_v2(df.query(f'{col}==True'))[s_ini:s_fin, 1, idx_cama]
            curva_no_elegible = curvaHosp_v2(df.query(f'{col}==False'))[s_ini:s_fin, 1, idx_cama]

            curvas_mb[cama]['elegible'].append((k, year, curva_elegible))
            curvas_mb[cama]['no_elegible'].append((k, year, curva_no_elegible))
df_lrti_vrs_MB = pd.concat([
    construir_df(curvas_mb, elegible=True, cama='MB'),
    construir_df(curvas_mb, elegible=False, cama='MB')
])
df_lrti_vrs_MB.to_csv(f'{folder_path}/Paper/FINALES/df_lrti_vrs_MB_gonzalo.csv')


############################ Egresos LRTI Amal ############################
egresos_lrti_amal = diseases_amal[(diseases_amal['irag'])|(diseases_amal['vrs1'])]
for año in [2019, 2022, 2023]:
    col_catchup = f'catchup2{año}'
    col_posible = f'posible_{año}'
    col_inseason = f'inseason_{año}'
    egresos_lrti_amal[col_catchup] = egresos_lrti_amal.apply(
        lambda row: df.catchup(row['fechaNac'], año) if row[col_posible] and not row[col_inseason] else -1,
        axis=1
    )
for y in [2019,2022,2023]:
    for k in range(0,13):
        col = f'elegible{k}_{y}'
        col_catchup = f'catchup2{y}'
        col_posible = f'posible_{y}'
        col_inseason = f'inseason_{y}'
        col_riesgo = f'riesgo_{y}'
        egresos_lrti_amal[col] =  ((egresos_lrti_amal[col_riesgo]) | (egresos_lrti_amal[col_inseason]) |
                               ((egresos_lrti_amal[col_catchup]>=0) & (egresos_lrti_amal[col_catchup]<=k)))
                
egresos_lrti_estadias_amal = transform_translations(egresos_lrti_amal,cf.areasUPC,cf.areasMB,min_fecha)
egresos_lrti_estadias_2019_amal = egresos_lrti_estadias_amal.query('indexIng>=@season[0][0]').query('indexIng<=@season[0][1]')
egresos_lrti_estadias_2022_amal = egresos_lrti_estadias_amal.query('indexIng>=@season[3][0]').query('indexIng<=@season[3][1]')
egresos_lrti_estadias_2023_amal = egresos_lrti_estadias_amal.query('indexIng>=@season[4][0]').query('indexIng<=@season[4][1]')


años = [2019, 2022, 2023]
campos = {'UPC': 1, 'MB': 0}
rango_k = range(13)
season_map = {2019: 0, 2022: 3, 2023: 4}
egresos_map = {
    2019: egresos_lrti_estadias_2019_amal,
    2022: egresos_lrti_estadias_2022_amal,
    2023: egresos_lrti_estadias_2023_amal,
}
### UPC VRS
curvas_upc = {cama: {'elegible': [], 'no_elegible': []} for cama in campos}
for cama, idx_cama in campos.items():
    for year in años:
        df = egresos_map[year]
        s_ini, s_fin = season[season_map[year]]

        for k in rango_k:
            col = f'elegible{k}_{year}'

            curva_elegible = curvaHosp_v2(df.query(f'{col}==True'))[s_ini:s_fin, 1, idx_cama]
            curva_no_elegible = curvaHosp_v2(df.query(f'{col}==False'))[s_ini:s_fin, 1, idx_cama]

            curvas_upc[cama]['elegible'].append((k, year, curva_elegible))
            curvas_upc[cama]['no_elegible'].append((k, year, curva_no_elegible))
df_lrti_vrs_UPC_amal = pd.concat([
    construir_df(curvas_upc, elegible=True, cama='UPC'),
    construir_df(curvas_upc, elegible=False, cama='UPC')
])
df_lrti_vrs_UPC_amal.to_csv(f'{folder_path}/Paper/FINALES/df_lrti_vrs_UPC_amal.csv')
### MB VRS
curvas_mb = {cama: {'elegible': [], 'no_elegible': []} for cama in campos}
for cama, idx_cama in campos.items():
    for year in años:
        df = egresos_map[year]
        s_ini, s_fin = season[season_map[year]]

        for k in rango_k:
            col = f'elegible{k}_{year}'

            curva_elegible = curvaHosp_v2(df.query(f'{col}==True'))[s_ini:s_fin, 1, idx_cama]
            curva_no_elegible = curvaHosp_v2(df.query(f'{col}==False'))[s_ini:s_fin, 1, idx_cama]

            curvas_mb[cama]['elegible'].append((k, year, curva_elegible))
            curvas_mb[cama]['no_elegible'].append((k, year, curva_no_elegible))
df_lrti_vrs_MB_amal = pd.concat([
    construir_df(curvas_mb, elegible=True, cama='MB'),
    construir_df(curvas_mb, elegible=False, cama='MB')
])
df_lrti_vrs_MB_amal.to_csv(f'{folder_path}/Paper/FINALES/df_lrti_vrs_MB_amal.csv')


###############################################################################
############################### Atenciones de Urgencias #######################
###############################################################################
# Se calcula la proporción de público y privado con la información de los hospitales de vigilancia
virus = (pd.read_csv(f'{folder_path_input}/IMPACTO/virus-edad.csv')
            .assign(total = lambda x: x['0a1']+x['1a4']+x['15a54']+x['5a14']+x['55a64']+x['65+'],
                            de5 = lambda x: x['15a54']+x['5a14']+x['55a64']+x['65+'])
                    .replace(['Influenza A','Influenza B'],'Influenza')
                    .groupby(['year','week','Agente'])
                    .agg({'0a1':'sum','1a4':'sum','de5':'sum'})
                    .reset_index()
                    .melt(id_vars=['year','week','Agente'], 
                                value_vars=['0a1','1a4','de5'],
                                var_name='age_group', 
                                value_name='count')
                    .pivot_table(index=['year','week','age_group'], columns='Agente', values='count', fill_value=0)
                    .reset_index())
virus_01 = (virus[virus['age_group'].isin(['0a1'])]
         .groupby(['year','week'])
         .agg({
             'Adenovirus':'sum','Influenza':'sum','Metapneumovirus':'sum','Parainfluenza':'sum','Virus Respiratorio Sincicial':'sum'
         })
         .reset_index()
         .assign(virus_total = lambda x: x['Adenovirus']+x['Influenza']+x['Metapneumovirus']+x['Parainfluenza']+x['Virus Respiratorio Sincicial'])
         .assign(vrs_prop_01 = lambda x: x['Virus Respiratorio Sincicial']/x['virus_total'])
         .rename(columns={'week':'epiweek'})
         .fillna(0))
virus_14 = (virus[virus['age_group'].isin(['1a4'])]
         .groupby(['year','week'])
         .agg({
             'Adenovirus':'sum','Influenza':'sum','Metapneumovirus':'sum','Parainfluenza':'sum','Virus Respiratorio Sincicial':'sum'
         })
         .reset_index()
         .assign(virus_total = lambda x: x['Adenovirus']+x['Influenza']+x['Metapneumovirus']+x['Parainfluenza']+x['Virus Respiratorio Sincicial'])
         .assign(vrs_prop_14 = lambda x: x['Virus Respiratorio Sincicial']/x['virus_total'])
         .rename(columns={'week':'epiweek'})
         .fillna(0))
atenciones = (pd
              .read_parquet(f'{folder_path}/Paper/IMPACTO/urgencias.parquet')
              .query('~categoria_ICD10.isin(["Upper respiratory tract infection"])')
              .assign(fechaIng = lambda df: pd.to_datetime({'year': df['year'], 'month': df['month'], 'day': df['day']}, format='%Y-%m-%d'))
              .assign(epiweek = lambda x: x['fechaIng'].dt.isocalendar().week))
atenciones2019 = atenciones.query('fechaIng>=@date_inicio_2019').query('fechaIng<=@date_fin_2019')
atenciones2022 = atenciones.query('fechaIng>=@date_inicio_2022').query('fechaIng<=@date_fin_2022')
atenciones2023 = atenciones.query('fechaIng>=@date_inicio_2023').query('fechaIng<=@date_fin_2023')
atenciones_lrti = (pd
                   .concat([atenciones2019, atenciones2022, atenciones2023], ignore_index=True)
                   .merge(virus_01[['year','epiweek','vrs_prop_01']],on=['year','epiweek'],how='left'))
egresos_g = pd.read_parquet(f'{folder_path}/Paper/IMPACTO/egresos_gonzalo.parquet').query('irag==True')
egresos_g = pd.concat([egresos_g.query('fechaIng>=@date_inicio_2019').query('fechaIng<=@date_fin_2019'),
                       egresos_g.query('fechaIng>=@date_inicio_2022').query('fechaIng<=@date_fin_2022'),
                       egresos_g.query('fechaIng>=@date_inicio_2023').query('fechaIng<=@date_fin_2023')])
egresos_a = pd.read_parquet(f'{folder_path}/Paper/IMPACTO/egresos_amal_04.parquet').query('irag==True')
egresos_a = pd.concat([egresos_a.query('fechaIng>=@date_inicio_2019').query('fechaIng<=@date_fin_2019'),
                       egresos_a.query('fechaIng>=@date_inicio_2022').query('fechaIng<=@date_fin_2022'),
                       egresos_a.query('fechaIng>=@date_inicio_2023').query('fechaIng<=@date_fin_2023')])

# Aquí calculamos distinto la proporción entre público y privado 
vrs01_g = (egresos_g[egresos_g['pneumonia']|egresos_g['acute_bronchitis_bronchiolitis']]
            .query('ageGroup=="menor_1"')
            .groupby(['ANO_ING','epiweek','vrs1'])
            .agg({'RUN':'count'})
            .reset_index())
vrs01_total_g = (egresos_g[egresos_g['pneumonia']|egresos_g['acute_bronchitis_bronchiolitis']]
                    .query('ageGroup=="menor_1"')
                    .groupby(['ANO_ING','epiweek'])
                    .agg({'RUN':'count'})
                    .reset_index()
                    .rename(columns={'RUN':'total'}))
vrs01_g = ((vrs01_g
           .query('vrs1==True')
           .merge(vrs01_total_g,how='left',on=['ANO_ING','epiweek'])
           .assign(vrs_prop = lambda x: x['RUN']/x['total']))[['ANO_ING','epiweek','vrs_prop']]
           .rename(columns={'ANO_ING':'year','epiweek':'week','vrs_prop':'vrs_prop_g'}))

vrs01_a = (egresos_a[egresos_a['pneumonia']|egresos_a['acute_bronchitis_bronchiolitis']]
            .query('ageGroup=="menor_1"')
            .groupby(['ANO_ING','epiweek','vrs1'])
            .agg({'RUN':'count'})
            .reset_index())
vrs01_total_a = (egresos_a[egresos_a['pneumonia']|egresos_a['acute_bronchitis_bronchiolitis']]
                    .query('ageGroup=="menor_1"')
                    .groupby(['ANO_ING','epiweek'])
                    .agg({'RUN':'count'})
                    .reset_index()
                    .rename(columns={'RUN':'total'}))
vrs01_a = ((vrs01_a
           .query('vrs1==True')
           .merge(vrs01_total_a,how='left',on=['ANO_ING','epiweek'])
           .assign(vrs_prop = lambda x: x['RUN']/x['total']))[['ANO_ING','epiweek','vrs_prop']]
           .rename(columns={'ANO_ING':'year','epiweek':'week','vrs_prop':'vrs_prop_a'}))

atenciones_lrti_total = (atenciones_lrti
                         .groupby(['year','epiweek'])
                         .agg({'Menores_1':'sum','De_1_a_4':'sum'})
                         .reset_index()
                         .rename(columns={'epiweek':'week'})
                         .fillna(0))
atenciones_lrti_vrs = (atenciones_lrti
                       .query('categoria_ICD10.isin(["Acute Bronchitis/Bronchiolitis","Pneumonia"])')
                       .assign(vrs_isp = lambda x: np.round(x['Menores_1']*x['vrs_prop_01']))
                       .groupby(['year','epiweek'])
                       .agg({'Menores_1':'sum','De_1_a_4':'sum','vrs_isp':'sum'})
                       .reset_index()
                       .rename(columns={'epiweek':'week'})
                       .fillna(0)
                       .merge(vrs01_a,how='left',on=['year','week'])
                       .merge(vrs01_g,how='left',on=['year','week'])
                       .assign(vrs_isp_a = lambda x: np.round(x['Menores_1']*x['vrs_prop_a']))
                       .assign(vrs_isp_g = lambda x: np.round(x['Menores_1']*x['vrs_prop_g']))
                       .drop(columns=['vrs_prop_a','vrs_prop_g']))

# ---------------------
# Versión gonzalo
# ---------------------
nac_ped_lambda = egresos_g[
    (egresos_g['ANO_ING'] >= 2022) &
    (egresos_g['vrs1']) &
    (egresos_g['edad_meses'] <= 48) &
    (egresos_g['edad_meses'] > 12) 
]
nac_ped_lambda = nac_ped_lambda.assign(
    edad = nac_ped_lambda['edad_meses'].apply(lambda x: "12a21" if x <= 21 else "21a48")
)
lambda_table = (nac_ped_lambda
    .groupby('edad')
    .size()
    .reset_index(name='n')
)
lambda_table['prop'] = lambda_table['n'] / lambda_table['n'].sum()
lambda_value = lambda_table[lambda_table['edad'] == "12a21"]['prop'].iloc[0]
atenciones_lrti_vrs = (
    atenciones_lrti_vrs
    .assign(vrs_total_g = lambda x: x['vrs_isp_g']+np.round(x['De_1_a_4']*lambda_value)))

# ---------------------
# Versión Amal
# ---------------------
# lambda 2019
# ---------------------
nac_ped_lambda = egresos_a[
    (egresos_a['ANO_ING'] == 2019) &
    (egresos_a['vrs1']) &
    (egresos_a['edad_meses'] <= 48) &
    (egresos_a['edad_meses'] > 12) 
]

nac_ped_lambda = nac_ped_lambda.assign(
    edad = nac_ped_lambda['edad_meses'].apply(lambda x: "12a13" if x <= 13 
                                              else ("13a14" if x <= 14
                                                    else ("14a15" if x <= 15
                                                           else ("15a16" if x <= 16
                                                                 else ("16a17" if x <= 17
                                                                       else ("17a18" if x <= 18
                                                                             else "18+")))))))
lambda_table = (nac_ped_lambda
    .groupby('edad')
    .size()
    .reset_index(name='n')
)
lambda_table['prop'] = lambda_table['n'] / lambda_table['n'].sum()

lambdas_amal = []
semanas_amal = []
for i in range(13,41):
    semanas_amal.append(i)
    if i>=13 and i<=16:
        lambda_semana = lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
    elif i>=17 and i<=21:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0])
    elif i>=22 and i<=25:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0])
    elif i>=26 and i<=30:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "15a16"]['prop'].iloc[0])
    elif i>=31 and i<=35:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "15a16"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "16a17"]['prop'].iloc[0])
    elif i>=36 and i<=40:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "15a16"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "16a17"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "17a18"]['prop'].iloc[0])
    lambdas_amal.append(lambda_semana)

lambdas_amal_final_2019 = pd.DataFrame({'week':semanas_amal,'lambda':lambdas_amal})
lambdas_amal_final_2019['year']=2019
# ---------------------
# lambda 2022
# ---------------------
nac_ped_lambda = egresos_a[
    (egresos_a['ANO_ING'] == 2022) &
    (egresos_a['vrs1']) &
    (egresos_a['edad_meses'] <= 48) &
    (egresos_a['edad_meses'] > 12) 
]

nac_ped_lambda = nac_ped_lambda.assign(
    edad = nac_ped_lambda['edad_meses'].apply(lambda x: "12a13" if x <= 13 
                                              else ("13a14" if x <= 14
                                                    else ("14a15" if x <= 15
                                                           else ("15a16" if x <= 16
                                                                 else ("16a17" if x <= 17
                                                                       else ("17a18" if x <= 18
                                                                             else "18+")))))))
lambda_table = (nac_ped_lambda
    .groupby('edad')
    .size()
    .reset_index(name='n')
)
lambda_table['prop'] = lambda_table['n'] / lambda_table['n'].sum()

lambdas_amal = []
semanas_amal = []
for i in range(13,41):
    semanas_amal.append(i)
    if i>=13 and i<=16:
        lambda_semana = lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
    elif i>=17 and i<=21:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0])
    elif i>=22 and i<=25:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0])
    elif i>=26 and i<=30:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "15a16"]['prop'].iloc[0])
    elif i>=31 and i<=35:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "15a16"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "16a17"]['prop'].iloc[0])
    elif i>=36 and i<=40:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "15a16"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "16a17"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "17a18"]['prop'].iloc[0])
    lambdas_amal.append(lambda_semana)

lambdas_amal_final_2022 = pd.DataFrame({'week':semanas_amal,'lambda':lambdas_amal})
lambdas_amal_final_2022['year']=2022
# ---------------------
# lambda 2023
# ---------------------
nac_ped_lambda = egresos_a[
    (egresos_a['ANO_ING'] >= 2022) &
    (egresos_a['vrs1']) &
    (egresos_a['edad_meses'] <= 48) &
    (egresos_a['edad_meses'] > 12) 
]

nac_ped_lambda = nac_ped_lambda.assign(
    edad = nac_ped_lambda['edad_meses'].apply(lambda x: "12a13" if x <= 13 
                                              else ("13a14" if x <= 14
                                                    else ("14a15" if x <= 15
                                                           else ("15a16" if x <= 16
                                                                 else ("16a17" if x <= 17
                                                                       else ("17a18" if x <= 18
                                                                             else "18+")))))))
lambda_table = (nac_ped_lambda
    .groupby('edad')
    .size()
    .reset_index(name='n')
)
lambda_table['prop'] = lambda_table['n'] / lambda_table['n'].sum()

lambdas_amal = []
semanas_amal = []
for i in range(13,41):
    semanas_amal.append(i)
    if i>=13 and i<=16:
        lambda_semana = lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
    elif i>=17 and i<=21:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0])
    elif i>=22 and i<=25:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0])
    elif i>=26 and i<=30:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "15a16"]['prop'].iloc[0])
    elif i>=31 and i<=35:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "15a16"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "16a17"]['prop'].iloc[0])
    elif i>=36 and i<=40:
        lambda_semana = (lambda_table[lambda_table['edad'] == "12a13"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "13a14"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "14a15"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "15a16"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "16a17"]['prop'].iloc[0]
                         +lambda_table[lambda_table['edad'] == "17a18"]['prop'].iloc[0])
    lambdas_amal.append(lambda_semana)

lambdas_amal_final_2023 = pd.DataFrame({'week':semanas_amal,'lambda':lambdas_amal})
lambdas_amal_final_2023['year']=2023

lambdas_amal_final = pd.concat([lambdas_amal_final_2019,lambdas_amal_final_2022,lambdas_amal_final_2023],ignore_index=True)
atenciones_lrti_vrs_final = (atenciones_lrti_vrs
                                .merge(lambdas_amal_final,how='left',on=['year','week'])
                                .assign(vrs_total_a = lambda x: x['vrs_isp_a']+np.round(x['De_1_a_4']*x['lambda']))
                                )[['year','week','vrs_total_g','vrs_total_a']]

############################ Curvas finales de Atenciones de Urgencia ############################
# Originales públicas
atenciones_lrti_vrs_final.to_csv(f'{folder_path}/Paper/FINALES/atenciones_publicas_gonzalo_y_amal.csv')
# Públicas y privadas
prop_pub_priv_a = (egresos_a
                        [(egresos_a['ANO_ING'].isin([2019,2022,2023])) &
                        (egresos_a['vrs1']) &
                        (egresos_a['edad_meses'] <= 48) ]
                        .groupby(['ANO_ING','origen_estab'])
                        .agg({'RUN':'count'})
                        .reset_index()
                        .rename(columns={'RUN':'n'})
                        .assign(prop_pub_a = lambda x: x['n']/x.groupby('ANO_ING')['n'].transform('sum'))
                        .drop(columns=['n']))
prop_pub_priv_g = (egresos_g
                        [(egresos_g['ANO_ING'].isin([2019,2022,2023])) &
                        (egresos_g['vrs1']) &
                        (egresos_g['edad_meses'] <= 48) ]
                        .groupby(['ANO_ING','origen_estab'])
                        .agg({'RUN':'count'})
                        .reset_index()
                        .rename(columns={'RUN':'n'})
                        .assign(prop_pub_g = lambda x: x['n']/x.groupby('ANO_ING')['n'].transform('sum'))
                        .drop(columns=['n']))
prop_pub = (prop_pub_priv_a
            .merge(prop_pub_priv_g,how='left',on=['ANO_ING','origen_estab'])
            .rename(columns={'ANO_ING':'year'})
            .query('origen_estab=="public"')
            .drop(columns={'origen_estab'}))
atenciones_lrti_vrs_final_pub_priv = (atenciones_lrti_vrs_final
                                        .merge(prop_pub,how='left',on=['year'])
                                        .assign(vrs_pub_priv_g = lambda x: np.round(x['vrs_total_g']/x['prop_pub_g']))
                                        .assign(vrs_pub_priv_a = lambda x: np.round(x['vrs_total_a']/x['prop_pub_a']))
                                        )
atenciones_lrti_vrs_final_pub_priv.to_csv(f'{folder_path}/Paper/FINALES/atenciones_pub_y_priv_gonzalo_y_amal.csv')

###############################################################################
###############################################################################
###############################################################################
################################# Simulación  #################################
###############################################################################
###############################################################################
###############################################################################

#################### Generación de matriz MA-HA-ICU ####################
upc_a = (pd.read_csv(f'{folder_path}/Paper/FINALES/df_lrti_vrs_UPC_amal.csv')).drop(columns=['Unnamed: 0'])
upc_g = (pd.read_csv(f'{folder_path}/Paper/FINALES/df_lrti_vrs_UPC_gonzalo.csv')).drop(columns=['Unnamed: 0'])
mb_a = (pd.read_csv(f'{folder_path}/Paper/FINALES/df_lrti_vrs_MB_amal.csv')).drop(columns=['Unnamed: 0'])
mb_g = (pd.read_csv(f'{folder_path}/Paper/FINALES/df_lrti_vrs_MB_gonzalo.csv')).drop(columns=['Unnamed: 0'])
at = (pd.read_csv(f'{folder_path}/Paper/FINALES/atenciones_pub_y_priv_gonzalo_y_amal.csv')).drop(columns=['Unnamed: 0'])
curva_nacimientos = pd.read_excel(f'{folder_path}/Paper/FINALES/Tablas/curva_nacimientos.xlsx')
columnas_dia = [str(i) for i in range(182)]  
# ---------------------------
# ATENCIONES URGENCIA MA
# ---------------------------
# Atenciones urgencias gonzalo
at_g = (at
        .groupby(['year'])
        .agg({'vrs_pub_priv_g':'sum'})
        .reset_index()
        .rename(columns={'vrs_pub_priv_g':'MA'}))
# Atenciones urgencias amal
at_a = (at
        .groupby(['year'])
        .agg({'vrs_pub_priv_a':'sum'})
        .reset_index()
        .rename(columns={'vrs_pub_priv_a':'MA'}))
# Aproximación Atenciones por Estrategías 
# por eso yo calcule más lambdas, porque según la estrategía depende el valor por el cual se debe extrapolar 
# de las atenciones de urgencias públicas al total
for i in range(13):
    col_estrategia = f'estrategia_{i}'
    col_catchup = f'catchup{i}'
    curva_nacimientos[col_estrategia] = (curva_nacimientos['inseason']+curva_nacimientos[col_catchup])/curva_nacimientos['posibles']
estrategia = []
proporcion = []
years = []
for y in [2019,2022,2023]:
    for i in range(13):
        col_estrategia = f'estrategia_{i}'
        proporcion.append(curva_nacimientos[curva_nacimientos['year']==y][col_estrategia].values[0])
        estrategia.append(i)
        years.append(y)

prop_est_year = pd.DataFrame({'year':years,'estrategia':estrategia,'proporcion':proporcion})
aten_g = (prop_est_year
          .merge(at_g, on='year', how='left')
          .assign(MA = lambda x: np.round(x['MA']*x['proporcion']))
          .drop(columns=['proporcion']))
aten_a = (prop_est_year
          .merge(at_a, on='year', how='left')
          .assign(MA = lambda x: np.round(x['MA']*x['proporcion']))
          .drop(columns=['proporcion']))

# ---------------------------
# UPC (ICU) y MB (HA)
# ---------------------------
# Datos Formato Gonzalo
upc_largo = upc_g.melt(
    id_vars=["estrategia", "year",'elegible'],     
    value_vars=columnas_dia,           
    var_name="dia",                    
    value_name="ICU"                  
)
upc_largo["dia"] = upc_largo["dia"].astype(int)
mb_largo = mb_g.melt(
    id_vars=["estrategia", "year",'elegible'],     
    value_vars=columnas_dia,           
    var_name="dia",                    
    value_name="HA"                  
)
mb_largo["dia"] = mb_largo["dia"].astype(int)
total_g = upc_largo.merge(mb_largo,on=['dia','estrategia','year','elegible'],how='left')

total_g = (total_g
           .query('elegible == 1')
           .drop(columns=['elegible'])
           .groupby(['estrategia','year'])
           .agg({'ICU':'sum','HA':'sum'})
           .reset_index())
total_g = (total_g
            .merge(aten_g,on=['year','estrategia'],how='left'))
total_g.to_csv(f'{folder_path}/Paper/FINALES/MA-HA-ICU_g.csv')
# Datos Formato Amal 
upc_largo = upc_a.melt(
    id_vars=["estrategia", "year",'elegible'],     
    value_vars=columnas_dia,           
    var_name="dia",                    
    value_name="ICU"                  
)
upc_largo["dia"] = upc_largo["dia"].astype(int)
mb_largo = mb_a.melt(
    id_vars=["estrategia", "year",'elegible'],     
    value_vars=columnas_dia,           
    var_name="dia",                    
    value_name="HA"                  
)
mb_largo["dia"] = mb_largo["dia"].astype(int)
total_a= upc_largo.merge(mb_largo,on=['dia','estrategia','year','elegible'],how='left')
total_a = (total_a
           .groupby(['estrategia','year'])
           .agg({'ICU':'sum','HA':'sum'})
           .reset_index())
total_a = (total_a
            .merge(aten_a,on=['year','estrategia'],how='left'))
total_a.to_csv(f'{folder_path}/Paper/FINALES/MA-HA-ICU_a.csv')

#################### SIMULACION ####################
total_g = pd.read_csv(f'{folder_path}/Paper/FINALES/MA-HA-ICU_g.csv').drop(columns=['Unnamed: 0'])
total_a = pd.read_csv(f'{folder_path}/Paper/FINALES/MA-HA-ICU_a.csv').drop(columns=['Unnamed: 0'])
np.random.seed(42)
# alpha calculado para replicar valores de eficiencia
# ma_mean = 76.4
# ma_min = 62.3
# ma_max = 85.2
# ha_mean = 76.8
# ha_min = 49.4
# ha_max = 89.4
# hau_mean = 78.6
# hau_min = 48.8
# hau_max = 91
np.random.seed(42)
alpha = [39.60, 0.21, 0.93, 11.09]  
def simulate_prevention(row, eff_MA, eff_HA, eff_ICU):
    new_row = row.copy()

    # MA: prevenir con probabilidad eff_MA
    prevented_MA = np.random.binomial(row['MA'], eff_MA)
    # MA: restante MA es MA menos lo prevenido
    new_row['MA'] -= prevented_MA

    # HA: prevenir con probabilidad eff_MA
    prevented_HA = np.random.binomial(row['HA'], eff_MA)
    # HA: restante HA es HA menos lo prevenido
    remaining_HA = row['HA'] - prevented_HA
    # HA: transformado a MA
    converted_HA_to_MA = np.random.binomial(remaining_HA, (eff_HA - eff_MA)/(1 - eff_MA))
    # HA: restante HA es restante HA menos lo transformado a MA
    new_row['HA'] = remaining_HA - converted_HA_to_MA
    # HA: se agrega a MA lo transformado de HA a MA
    new_row['MA'] += converted_HA_to_MA

    # ICU: prevenir con probabilidad eff_MA
    prevented_ICU = np.random.binomial(row['ICU'], eff_MA)
    # ICU: restante ICU es ICU menos lo prevenido
    remaining_ICU = row['ICU'] - prevented_ICU
    # ICU: transformado a MA
    converted_ICU_to_MA = np.random.binomial(remaining_ICU, (eff_HA - eff_MA)/(1 - eff_MA))
    # ICU: restante ICU es ICU menos lo transformado a MA
    remaining_ICU -= converted_ICU_to_MA
    # ICU: se agrega a MA lo transformado de ICU a MA
    new_row['MA'] += converted_ICU_to_MA
    # ICU: restante ICU transformado a HA
    converted_ICU_to_HA = np.random.binomial(remaining_ICU, (eff_ICU - eff_HA)/(1 - eff_HA))
    # ICU: restante ICU es restance ICU menos lo transformado a HA
    remaining_ICU -= converted_ICU_to_HA
    # ICU: rse agrega a HA lo transformado de ICU a HA
    new_row['HA'] += converted_ICU_to_HA
    # ICU final
    new_row['ICU'] = remaining_ICU

    return new_row


# Simulacion gonzalo
N = 1000
e = np.random.dirichlet(alpha,N)
resultados = []
for i in range(N):
    print(f"Iteración {i}/{N}")
    sim_result = total_g.query('elegible==1').apply(simulate_prevention, axis=1, args=(e[i][0],e[i][0]+e[i][1],e[i][0]+e[i][1]+e[i][2])).assign(sim_id=i)
    resultados.append(sim_result[['year','dia','estrategia', 'MA', 'HA', 'ICU', 'sim_id']])
df_montecarlo = pd.concat(resultados[1000:], ignore_index = True)
df_montecarlo.to_parquet(f'{folder_path}/Paper/FINALES/df_montecarlo_g.csv', engine='pyarrow', compression='snappy')
col='mean'
mean_mc = (df_montecarlo
        .groupby(['year','dia','estrategia'])
        .agg({'MA':col,'HA':col,'ICU':col})
        .reset_index()
        .rename(columns={'MA':f'MA_{col}','HA':f'HA_{col}','ICU':f'ICU{col}'}))
col='std'
std_mc = (df_montecarlo
        .groupby(['year','dia','estrategia'])
        .agg({'MA':col,'HA':col,'ICU':col})
        .reset_index()
        .rename(columns={'MA':f'MA_{col}','HA':f'HA_{col}','ICU':f'ICU{col}'}))
col='min'
min_mc = (df_montecarlo
        .groupby(['year','dia','estrategia'])
        .agg({'MA':col,'HA':col,'ICU':col})
        .reset_index()
        .rename(columns={'MA':f'MA_{col}','HA':f'HA_{col}','ICU':f'ICU{col}'}))
col='max'
max_mc = (df_montecarlo
        .groupby(['year','dia','estrategia'])
        .agg({'MA':col,'HA':col,'ICU':col})
        .reset_index()
        .rename(columns={'MA':f'MA_{col}','HA':f'HA_{col}','ICU':f'ICU{col}'}))
resumen = (mean_mc
           .merge(std_mc,on=['year','dia','estrategia'],how='left')
           .merge(min_mc,on=['year','dia','estrategia'],how='left')
           .merge(max_mc,on=['year','dia','estrategia'],how='left'))
resumen.to_csv(f'{folder_path}/Paper/FINALES/resumen_montecarlo_g.csv')

# Simulacion amal
N = 1000
e = np.random.dirichlet(alpha,N)
resultados = []
for i in range(N):
    print(f"Iteración {i}/{N}")
    sim_result = total_g.query('elegible==1').apply(simulate_prevention, axis=1, args=(e[i][0],e[i][0]+e[i][1],e[i][0]+e[i][1]+e[i][2])).assign(sim_id=i)
    resultados.append(sim_result[['year','dia','estrategia', 'MA', 'HA', 'ICU', 'sim_id']])
df_montecarlo = pd.concat(resultados[1000:], ignore_index = True)
df_montecarlo.to_parquet(f'{folder_path}/Paper/FINALES/df_montecarlo_a.csv', engine='pyarrow', compression='snappy')
col='mean'
mean_mc = (df_montecarlo
        .groupby(['year','dia','estrategia'])
        .agg({'MA':col,'HA':col,'ICU':col})
        .reset_index()
        .rename(columns={'MA':f'MA_{col}','HA':f'HA_{col}','ICU':f'ICU{col}'}))
col='std'
std_mc = (df_montecarlo
        .groupby(['year','dia','estrategia'])
        .agg({'MA':col,'HA':col,'ICU':col})
        .reset_index()
        .rename(columns={'MA':f'MA_{col}','HA':f'HA_{col}','ICU':f'ICU{col}'}))
col='min'
min_mc = (df_montecarlo
        .groupby(['year','dia','estrategia'])
        .agg({'MA':col,'HA':col,'ICU':col})
        .reset_index()
        .rename(columns={'MA':f'MA_{col}','HA':f'HA_{col}','ICU':f'ICU{col}'}))
col='max'
max_mc = (df_montecarlo
        .groupby(['year','dia','estrategia'])
        .agg({'MA':col,'HA':col,'ICU':col})
        .reset_index()
        .rename(columns={'MA':f'MA_{col}','HA':f'HA_{col}','ICU':f'ICU{col}'}))
resumen = (mean_mc
           .merge(std_mc,on=['year','dia','estrategia'],how='left')
           .merge(min_mc,on=['year','dia','estrategia'],how='left')
           .merge(max_mc,on=['year','dia','estrategia'],how='left'))
resumen.to_csv(f'{folder_path}/Paper/FINALES/resumen_montecarlo_a.csv')

###############################################################################
###############################################################################
###############################################################################
################################### Costos  ###################################
###############################################################################
###############################################################################
###############################################################################

###############################################################################
############################### MA-HA-ICU #######################
###############################################################################
resumen_a = pd.read_csv(f'{folder_path}/Paper/FINALES/resumen_montecarlo_a.csv').drop(columns=['Unnamed: 0']).astype({'estrategia':int})
resumen_g = pd.read_csv(f'{folder_path}/Paper/FINALES/resumen_montecarlo_g.csv').drop(columns=['Unnamed: 0']).astype({'estrategia':int})
original_a = pd.read_csv(f'{folder_path}/Paper/FINALES/MA-HA-ICU_a.csv').drop(columns=['Unnamed: 0']).astype({'estrategia':int})
original_g = pd.read_csv(f'{folder_path}/Paper/FINALES/MA-HA-ICU_g.csv').drop(columns=['Unnamed: 0']).astype({'estrategia':int})
p_basic_bed = 414
p_icu_bed = 1083
emergency = 192
# Datos Gonzalo
mc = (resumen_g
        .groupby(['year','estrategia'])
        .agg({'MA_mean':'sum','HA_mean':'sum','ICUmean':'sum'})
        .reset_index()
        .assign(MA_mc = lambda x: np.round(x['MA_mean'],0),
                HA_mc = lambda x: np.round(x['HA_mean'],0),
                ICU_mc = lambda x: np.round(x['ICUmean'],0))
        .astype({'estrategia':int,'MA_mc':int,'HA_mc':int,'ICU_mc':int}))
final = (original_g
         .merge(mc[['year','estrategia','MA_mc','HA_mc','ICU_mc']], on=['year','estrategia'], how='left'))
final = (final
         .assign(
             ICU_costos = lambda x: x['ICU'] * p_icu_bed,
             MA_costos = lambda x: x['MA'] * emergency,
             HA_costos = lambda x: x['HA'] * p_basic_bed,
             costo_total = lambda x: x['ICU_costos'] + x['MA_costos'] + x['HA_costos'],
             ICU_mc_costos = lambda x: x['ICU_mc'] * p_icu_bed,
             MA_mc_costos = lambda x: x['MA_mc'] * emergency,
             HA_mc_costos = lambda x: x['HA_mc'] * p_basic_bed,
             costo_mc_total = lambda x: x['ICU_mc_costos'] + x['MA_mc_costos'] + x['HA_mc_costos']))
final.to_excel(f'{folder_path}/Paper/FINALES/costos_MA_HA_ICU_g.xlsx', index=False)
# Datos Amal
mc = (resumen_a
        .groupby(['year','estrategia'])
        .agg({'MA_mean':'sum','HA_mean':'sum','ICUmean':'sum'})
        .reset_index()
        .assign(MA_mc = lambda x: np.round(x['MA_mean'],0),
                HA_mc = lambda x: np.round(x['HA_mean'],0),
                ICU_mc = lambda x: np.round(x['ICUmean'],0))
        .astype({'estrategia':int,'MA_mc':int,'HA_mc':int,'ICU_mc':int}))
final = (original_a
         .merge(mc[['year','estrategia','MA_mc','HA_mc','ICU_mc']], on=['year','estrategia'], how='left'))
final = (final
         .assign(
             ICU_costos = lambda x: x['ICU'] * p_icu_bed,
             MA_costos = lambda x: x['MA'] * emergency,
             HA_costos = lambda x: x['HA'] * p_basic_bed,
             costo_total = lambda x: x['ICU_costos'] + x['MA_costos'] + x['HA_costos'],
             ICU_mc_costos = lambda x: x['ICU_mc'] * p_icu_bed,
             MA_mc_costos = lambda x: x['MA_mc'] * emergency,
             HA_mc_costos = lambda x: x['HA_mc'] * p_basic_bed,
             costo_mc_total = lambda x: x['ICU_mc_costos'] + x['MA_mc_costos'] + x['HA_mc_costos']))
final.to_excel(f'{folder_path}/Paper/FINALES/costos_MA_HA_ICU_a.xlsx', index=False)

###############################################################################
################################# Nacimientos #################################
###############################################################################
### Curva de Costos de Nacimientos
curva_costos_adquisicion = pd.DataFrame({'year':[2019,2022,2023]})
curva_costos_adquisicion['costo'] = 'adquisicion'
curva_costos_administracion = pd.DataFrame({'year':[2019,2022,2023]})
curva_costos_administracion['costo'] = 'administracion'
curva_costos = pd.DataFrame({'year':[2019,2022,2023]})
curva_costos['costo'] = 'Total'
curva_ahorros = pd.DataFrame({'year':[2019,2022,2023]})
cant_nirse = []
for k in range(0,13):
    col = f'catchup{k}'
    cant = (curva_nacimientos['riesgo_inseason']
            +curva_nacimientos['riesgo_catchup']
            +curva_nacimientos['inseason']
            +curva_nacimientos[col])
    cant_nirse.append(cant)
cant_adm = []
for k in range(0,13):
    col = f'catchup{k}'
    adm = (curva_nacimientos['riesgo_catchup']
            +curva_nacimientos[col])
    cant_adm.append(adm)
costo_pali = [cf.palivizumab_19,cf.palivizumab_22,cf.palivizumab_23]
for k in range(0,13):
    col = f'Estrategía {k}'
    curva_costos_adquisicion[col] = cant_nirse[k]*cf.nirsevimab
    curva_costos_administracion[col] = cant_adm[k]*cf.dose_administration
    curva_costos[col] = cant_adm[k]*cf.dose_administration+cant_nirse[k]*cf.nirsevimab
    curva_ahorros[col] = cant_nirse[k]*cf.nirsevimab+cant_adm[k]*cf.dose_administration-costo_pali
curva_costos_adquisicion
curva_costos_administracion
curva_costos = pd.concat([curva_costos_administracion,curva_costos_adquisicion,curva_costos])
curva_costos.to_excel(f'{folder_path}/Paper/FINALES/Tablas/curva_costos_nacimientos.xlsx')
curva_ahorros.to_excel(f'{folder_path}/Paper/FINALES/Tablas/curva_ahorros_nacimientos.xlsx')
