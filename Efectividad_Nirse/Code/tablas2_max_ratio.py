import os 
import sys

import pandas as pd
import numpy as np

import pickle

from preproces_prod3 import *
from matching_case_control import *

import inv
from IPython.core.display import display

local_path = os.getcwd()
code_root = os.path.abspath(os.path.join(local_path, '..', 'Code'))

if code_root not in sys.path:
    sys.path.insert(0, code_root)
code_root = os.path.abspath(os.path.join(local_path, '../..', 'inv'))

if code_root not in sys.path:
    sys.path.insert(0, code_root)


#all_data,analyze_vrs_data, integer_programming_matching_gurobi,match_nn_max_dist_weigths,comparar_medias_test, charly_mip, charly_double_mip, mylogit, results_match, tabla_marcel, tabla_final,summary_eff, summary_eff_aux


display.max_output_lines = 500  # Adjust the number as needed
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

np.random.seed(42)

path_actual = Path.cwd()
path_data = path_actual.parent / 'Data'

df_vrs_match_case = pd.read_csv(path_data/'df_vrs_match_case_11_04_25.csv', index_col=0)  #df_vrs_match_case_s39_2012
df_upc_match_case = pd.read_csv(path_data/'df_upc_match_case_11_04_25.csv', index_col=0) #df_upc_match_case_s39_2012

icd10_codes_fr = [
    "N390",   # Infección del tracto urinario (sin síntomas respiratorios) INFECCION DE VÍAS URINARIAS, SITIO NO ESPECIFICADO
    "A090",     # Gastroenteritis aguda (sin síntomas respiratorios) OTRAS GASTROENTERITIS Y COLITIS NO ESPECIFICADAS DE ORIGEN INFECCIOSO
    "A099",     # Gastroenteritis aguda (sin síntomas respiratorios) GASTROENTERITIS Y COLITIS DE ORIGEN NO ESPECIFICADO
    "A09X",     # Gastroenteritis aguda (sin síntomas respiratorios) DIARREA Y GASTROENTERITIS DE PRESUNO ORIGEN INFECCIOSO
    "R100",  # Cólico infantil (sin fiebre ni síntomas respiratorios) ABDOMEN AGUDO
    "R101",  # Cólico infantil (sin fiebre ni síntomas respiratorios) DOLOR ABDOMINAL LOCALIZADO EN PARTE SUPERIOR
    "R102",  # Cólico infantil (sin fiebre ni síntomas respiratorios) DOLOR PELVICO Y PERINEAL
    "R103",  # Cólico infantil (sin fiebre ni síntomas respiratorios) DOLOR LOCALIZADO EN OTRAS PARTES INFERIORES DEL ABDOMEN
    "R104",  # Cólico infantil (sin fiebre ni síntomas respiratorios) OTROS DOLORES ABDOMINALES Y LOS NO ESPECIFICADOS
    "R11X",  # Cólico infantil (sin fiebre ni síntomas respiratorios) NAUSEA Y VOMITO
    "R634",   # Pérdida de peso (sin fiebre ni síntomas respiratorios) PERDIDA ANORMAL DE PESO
    "R633",   # Dificultades de alimentación (sin fiebre ni síntomas respiratorios) DIFICULTADES Y MALA ADMINISTRACION DE LA ALIMENTACION
    "P599",   # Ictericia neonatal (sin fiebre ni síntomas respiratorios) ICTERICIA NEONATAL, NO ESPECIFICADA
    "R681",  # Llanto anormal (sin fiebre ni síntomas respiratorios) SINTOMAS NO ESPECIFICOS PROPIOS DE LA INFANCIA
    "S099",  # Traumatismo craneal (sin fiebre ni síntomas respiratorios)
    "Z539"    # Cirugía de emergencia (sin fiebre ni síntomas respiratorios)
    ################### AÑADIDOS POR MI #######################
    'A099', 
    'A080', 
    'A083', 
    'A082', 
    'A084',
    'A090', 
    'A081', 
    'A085'
]

ruts_cariopaticos_1 = pd.read_csv(path_data/'lista_pacientes_riesgo.csv').query('card1')
ruts_cariopaticos_2 = pd.read_csv(path_data/'lista_pacientes_riesgo.csv').query('card2')
ruts_displasia = pd.read_csv(path_data/'lista_pacientes_riesgo.csv').query('displ')

list_experiments=[]

df = (df_vrs_match_case
      .drop(columns=[col for col in df_vrs_match_case.columns if ((col.startswith('inm_')) | (col.startswith('MES_')) | (col.startswith('DIA_')) | 
                                               (col.startswith('SERC_')) | ((col.startswith('ANO_') & (col.endswith('TRAS')))) | 
                                               (col.startswith('AREAF')) | (col.startswith('dias_en')) | (col.startswith('fecha_tras')) |
                                               (col.endswith('Dall')) )])
      .assign(cardio1 = lambda x: x.RUN.isin(ruts_cariopaticos_1.RUN.unique()).astype(int),
              displa = lambda x: x.RUN.isin(ruts_displasia.RUN.unique()).astype(int),)
      .copy()
      )

df_upc = (df_upc_match_case
      .drop(columns=[col for col in df_upc_match_case.columns if ((col.startswith('inm_')) | (col.startswith('MES_')) | (col.startswith('DIA_')) | 
                                               (col.startswith('SERC_')) | ((col.startswith('ANO_') & (col.endswith('TRAS')))) | 
                                               (col.startswith('AREAF')) | (col.startswith('dias_en')) | (col.startswith('fecha_tras')) |
                                               (col.endswith('Dall')) )])
      .assign(cardio1 = lambda x: x.RUN.isin(ruts_cariopaticos_1.RUN.unique()).astype(int),
              displa = lambda x: x.RUN.isin(ruts_displasia.RUN.unique()).astype(int),)
      .copy()
      )

df.loc[df.RUN == 'ac483764636448d753930868dd3192a785e3728a2d81b3fc44dba45c3506255a', 'region'] = 'METROPOLITANA'
df.loc[df.COMUNA == 12202, 'region'] = 'MAGALLANES Y ANTARTICA'

#df_upc.loc[df.RUN == 'ac483764636448d753930868dd3192a785e3728a2d81b3fc44dba45c3506255a', 'region'] = 'METROPOLITANA'
#df_upc.loc[df.COMUNA == 12202, 'region'] = 'MAGALLANES Y ANTARTICA'

########################################################################################################################
match_vars_distance_nn=['edad_relativa','ingreso_relativo']
match_vars_exact_nn = ['group','region','PREVI','prematuro'] #NOMBRE_REGION
weights = {'edad_relativa': 1,'ingreso_relativo': 2} 

match_vars_distance_IP = ['edad_relativa']
match_vars_distance_IP_francia = ['edad_relativa','ingreso_relativo']

match_vars_exact_IP_previ = ['prematuro','region','group','PREVI'] # region
match_vars_exact_IP_not_previ = ['prematuro','region','group']


filtro_comba_dict = {  
                     'prematuros': '(prematuro == 1)',
                     'prematuros_seasonal': '(prematuro == 1)' + ' & ' + '(group == "SEASONAL")',
                     'prematuros_catchup': '(prematuro == 1)' + ' & ' + '(group == "CATCH_UP")',
                     
                    # 'VRS_general': 'MARCA==0',
                     
                     'displasia': 'RUN.isin(@ruts_displasia.RUN)',
                     'displasia_seasonal': '(RUN.isin(@ruts_displasia.RUN))' + ' & ' + '(group == "SEASONAL")',
                     'displasia_catchup': '(RUN.isin(@ruts_displasia.RUN))' + ' & ' + '(group == "CATCH_UP")',
                     
                     'cardiopats_1': 'RUN.isin(@ruts_cariopaticos_1.RUN)',
                     'cardiopats_1_seasonal': '(RUN.isin(@ruts_cariopaticos_1.RUN))' + ' & ' + '(group == "SEASONAL")',
                     'cardiopats_1_catchup': '(RUN.isin(@ruts_cariopaticos_1.RUN))' + ' & ' + '(group == "CATCH_UP")',
                    #  'cardiopats_2': 'RUN.isin(@ruts_cariopaticos_2.RUN)',

                      'cardiopats_O_prema': 'RUN.isin(@ruts_cariopaticos_1.RUN)' + ' | ' + '(prematuro == 1)',
                     'cardiopats_Y_prema': 'RUN.isin(@ruts_cariopaticos_1.RUN)' + ' & ' + '(prematuro == 1)',
                     
                    #  'cardiopats_Y_no_prema': 'RUN.isin(@ruts_cariopaticos_1.RUN)' + ' & ' + '(prematuro == 0)',
                    #  'prematuros_Y_no_cardio': '(prematuro == 1)' + ' & ' + '~RUN.isin(@ruts_cariopaticos_1.RUN)',
                     
                       'cardiopats_O_prematuros_catchup': '(RUN.isin(@ruts_cariopaticos_1.RUN)' + ' | ' + '(prematuro == 1))' + ' & ' + '(group == "CATCH_UP")',
                       'cardiopats_O_prematuros_seasonal': '(RUN.isin(@ruts_cariopaticos_1.RUN)' + ' | ' + '(prematuro == 1))' + ' & ' + '(group == "SEASONAL")',
                    #  'cardiopats_Y_prematuros_catchup': 'RUN.isin(@ruts_cariopaticos_1.RUN)' + ' & ' + '(prematuro == 1)' + ' & ' + '(group == "CATCH_UP")',
                    #  'cardiopats_Y_prematuros_seasonal': 'RUN.isin(@ruts_cariopaticos_1.RUN)' + ' & ' + '(prematuro == 1)' + ' & ' + '(group == "SEASONAL")',

                    }

grida_ratios=["1:5","1:7","1:10"] #["1:1","1:2","1:3","1:4","1:5","1:7","1:10"]

for i, ratiox in enumerate(grida_ratios):
    df_copy = df.copy()

    list_experiments_all_born = results_match(df_case_study=df_copy,
                                    df_control_study=df_copy,
                                    filtros_dic=filtro_comba_dict,
                                    match_vars_distance_nn=match_vars_distance_nn,
                                    match_vars_exact_nn=match_vars_exact_nn,
                                    match_vars_distance_IP=match_vars_distance_IP,
                                    match_vars_exact_IP=match_vars_exact_IP_not_previ,
                                    weights=weights, 
                                    list_experiments = [],
                                    nn=False,
                                    ratio=ratiox,
                                    covs = ['sexo','PESO','SEMANAS','cardio1','displa']) # ['sexo','SEMANAS','PESO','cardio1','prematuro']

    df_final = summary_eff_aux(list_experiments_all_born) #list_experiments_all_born,
    print(ratiox)
    print(df_final)
    # name_sheet = f"Tabla_2_ratio_{ratiox.split(':')[0]}a{ratiox.split(':')[1]}"
    # with pd.ExcelWriter(path_data/"resultados_comparados_copy.xlsx", engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    #     df_final.to_excel(writer, sheet_name=name_sheet)