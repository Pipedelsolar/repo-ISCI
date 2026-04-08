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

display.max_output_lines = 500  # Adjust the number as needed
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

np.random.seed(42)

path_actual = Path.cwd()
path_data = path_actual.parent / 'Data'

# df_vrs_match_case, _ = call_data('COHORTE_NIRSE_ACTUALIZADA_04_12_2025_ENCR.csv')  #   NAC_RNI_EGRESOS_ENTREGA_ISCI_11_04_2025_encr   COHORTE_NIRSE_ACTUALIZADA_04_12_2025_ENCR
df_vrs_match_case = pd.read_csv(path_data/'df_vrs_match_case_02_06_25_v2.csv') #df_vrs_match_case_02_06_25

# df_upc_match_case = pd.read_csv(path_data/'df_upc_match_case_11_04_25.csv', index_col=0) #df_upc_match_case_s39_2012

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

with open(path_data/'lista_ruts_cardio.pkl', 'rb') as f:
    lista_ruts_cardio = pickle.load(f)

with open(path_data/'lista_ruts_preterms.pkl', 'rb') as f:
    lista_ruts_preterms = pickle.load(f)
    
lista_risky = np.union1d(lista_ruts_cardio, lista_ruts_preterms)

# ruts_cariopaticos_1 = pd.read_csv(path_data/'lista_pacientes_riesgo.csv').query('card1')
# ruts_cariopaticos_2 = pd.read_csv(path_data/'lista_pacientes_riesgo.csv').query('card2')
ruts_displasia = pd.read_csv(path_data/'lista_pacientes_riesgo.csv').query('displ')

list_experiments=[]

df = (df_vrs_match_case
      .drop(columns=[col for col in df_vrs_match_case.columns if ((col.startswith('inm_')) | (col.startswith('MES_')) | (col.startswith('DIA_')) | 
                                               (col.startswith('SERC_')) | ((col.startswith('ANO_') & (col.endswith('TRAS')))) | 
                                               (col.startswith('AREAF')) | (col.startswith('dias_en')) | (col.startswith('fecha_tras')) |
                                               (col.endswith('Dall')) )])
      .assign(cardio1 = lambda x: x.RUN.isin(lista_ruts_cardio).astype(int),
              displa = lambda x: x.RUN.isin(ruts_displasia.RUN).astype(int),
              fecha_nac = lambda x: pd.to_datetime(x.fecha_nac, infer_datetime_format=True),
              year_nac = lambda x: x.fecha_nac.dt.year,
              month_nac = lambda x: x.fecha_nac.dt.month,
              criterio_pali_normal = lambda x: ((x.RUN.isin(lista_ruts_cardio)) | (x.SEMANAS<32) | (x.PESO<1500)).astype(int),
              pali = lambda x: np.where(x.criterio_pali_normal==1, 1, 
                                        np.where((x.year_nac==2023) & (x.month_nac.between(7, 9, inclusive='both')) & (x.SEMANAS<=34) & ((x.PESO<2500)), 1, 0)))
      .copy()
      )

# df_upc = (df_upc_match_case
#       .drop(columns=[col for col in df_upc_match_case.columns if ((col.startswith('inm_')) | (col.startswith('MES_')) | (col.startswith('DIA_')) | 
#                                                (col.startswith('SERC_')) | ((col.startswith('ANO_') & (col.endswith('TRAS')))) | 
#                                                (col.startswith('AREAF')) | (col.startswith('dias_en')) | (col.startswith('fecha_tras')) |
#                                                (col.endswith('Dall')) )])
#       .assign(cardio1 = lambda x: x.RUN.isin(lista_ruts_cardio).astype(int),
#               displa = lambda x: x.RUN.isin(ruts_displasia.RUN).astype(int),
#               fecha_nac = lambda x: pd.to_datetime(x.fecha_nac, infer_datetime_format=True),
#               year_nac = lambda x: x.fecha_nac.dt.year,
#               month_nac = lambda x: x.fecha_nac.dt.month,
#               criterio_pali_normal = lambda x: ((x.RUN.isin(lista_ruts_cardio)) | (x.SEMANAS<32) | (x.PESO<1500)).astype(int),
#               pali = lambda x: np.where(x.criterio_pali_normal==1, 1, 
#                                         np.where((x.year_nac==2023) & (x.month_nac.between(7, 9, inclusive='both')) & (x.SEMANAS<=34) & ((x.PESO<2500)), 1, 0)))
#       .copy()
#       )

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

match_vars_exact_IP_previ = ['prematuro','region','group','PREVI','cardio1'] # region
match_vars_exact_IP_not_previ = ['prematuro','region','group','cardio1']

dic_ratio = {  
                    #  'prematuros': '1:10',
                    #  'prematuros_seasonal': '1:1',
                    #  'prematuros_catchup': '1:10',
                     
                    # 'VRS_general': '1:7',
                    # 'UPC_general': '1:7',
                     
                    #  'displasia': 'RUN.isin(@ruts_displasia.RUN)',
                    #  'displasia_seasonal': '(RUN.isin(@ruts_displasia.RUN))' + ' & ' + '(group == "SEASONAL")',
                    #  'displasia_catchup': '(RUN.isin(@ruts_displasia.RUN))' + ' & ' + '(group == "CATCH_UP")',
                     
                    #  'cardiopats_1': '1:3',
                    #  'cardiopats_1_seasonal': '1:1',
                    #  'cardiopats_1_catchup': '1:3',
                     
                    #  'cardiopats_O_prema': '1:10', #1:2 igual resulta bien eff mas alto, pero mas ic
                    #  'cardiopats_Y_prema': '1:1',
                     
                    #  'cardiopats_O_prematuros_catchup': '1:10',
                    #  'cardiopats_O_prematuros_seasonal': '1:4',
                    }

########################################################## VRS ##########################################################

filtro_comba_dict_version_antigua = { 'High_risk': '((SEMANAS<32) | (PESO<1500) | (RUN.isin(@lista_ruts_cardio))) | ((32<=SEMANAS<=34) & (PESO>2500)) | ((PESO>=1500) & (SEMANAS==35)) & ~(RUN.isin(@lista_ruts_cardio))',
                      'High_risk_seasonal': '((SEMANAS<32) | (PESO<1500) | (RUN.isin(@lista_ruts_cardio))) | ((32<=SEMANAS<=34) & (PESO>2500)) | ((PESO>=1500) & (SEMANAS==35)) & ~(RUN.isin(@lista_ruts_cardio))' + ' & ' + '(group == "SEASONAL")',
                      'High_risk_catchup': '((SEMANAS<32) | (PESO<1500) | (RUN.isin(@lista_ruts_cardio))) | ((32<=SEMANAS<=34) & (PESO>2500)) | ((PESO>=1500) & (SEMANAS==35)) & ~(RUN.isin(@lista_ruts_cardio))' + ' & ' + '(group == "CATCH_UP")',  
                      
                      'pali': '(RUN.isin(@lista_ruts_cardio) | SEMANAS<32 | PESO<1500)',
                      'pali_seasonal': '(RUN.isin(@lista_ruts_cardio) | SEMANAS<32 | PESO<1500)' + ' & ' + '(group == "SEASONAL")',
                      'pali_catchup': '(RUN.isin(@lista_ruts_cardio) | SEMANAS<32 | PESO<1500)' + ' & ' + '(group == "CATCH_UP")',
                      
                      'super_prematuros': '(SEMANAS<32 | PESO<1500)',
                      'super_prematuros_seasonal': '(SEMANAS<32 | PESO<1500)' + ' & ' + '(group == "SEASONAL")', #DA NA
                      'super_prematuros_catchup': '(SEMANAS<32 | PESO<1500)' + ' & ' + '(group == "CATCH_UP")', #DA NA
                      
                    # 'VRS_general': 'MARCA==0',
                    
                    #  'displasia': 'RUN.isin(@ruts_displasia.RUN)',
                    #  'displasia_seasonal': '(RUN.isin(@ruts_displasia.RUN))' + ' & ' + '(group == "SEASONAL")',
                    #  'displasia_catchup': '(RUN.isin(@ruts_displasia.RUN))' + ' & ' + '(group == "CATCH_UP")',
                     
                      'cardiopats_1': 'RUN.isin(@lista_ruts_cardio)',
                      'cardiopats_1_seasonal': '(RUN.isin(@lista_ruts_cardio))' + ' & ' + '(group == "SEASONAL")',
                      'cardiopats_1_catchup': '(RUN.isin(@lista_ruts_cardio))' + ' & ' + '(group == "CATCH_UP")',
                    
                      'prematuros_no_pali' : '(SEMANAS<=35) & ~(RUN.isin(@lista_ruts_cardio) | (SEMANAS<32) | (PESO<1500))',
                      'prematuros_no_pali_seasonal' : '(SEMANAS<=35) & ~(RUN.isin(@lista_ruts_cardio) | SEMANAS<32 | PESO<1500)' + ' & ' + '(group == "SEASONAL")',
                      'prematuros_no_pali_catchup' : '(SEMANAS<=35) & ~(RUN.isin(@lista_ruts_cardio) | SEMANAS<32 | PESO<1500)' + ' & ' + '(group == "CATCH_UP")',
                      
                    #   'prematuros': '(SEMANAS<=35)', #'(prematuro == 1)',
                    #   'prematuros_seasonal': '(SEMANAS<=35)' + ' & ' + '(group == "SEASONAL")',
                    #   'prematuros_catchup': '(SEMANAS<=35)' + ' & ' + '(group == "CATCH_UP")',
                    #  'cardiopats_Y_prema': 'RUN.isin(@lista_ruts_cardio)' + ' & ' + '(prematuro == 1)',
                     
                    # #  'cardiopats_Y_no_prema': 'RUN.isin(@lista_ruts_cardio)' + ' & ' + '(prematuro == 0)',
                    # #  'prematuros_Y_no_cardio': '(prematuro == 1)' + ' & ' + '~RUN.isin(@lista_ruts_cardio)',
                     
                    #    'cardiopats_O_prematuros_catchup': '(RUN.isin(@lista_ruts_cardio)' + ' | ' + '(prematuro == 1))' + ' & ' + '(group == "CATCH_UP")',
                    #    'cardiopats_O_prematuros_seasonal': '(RUN.isin(@lista_ruts_cardio)' + ' | ' + '(prematuro == 1))' + ' & ' + '(group == "SEASONAL")',

                    }

filtro_comba_dict_labels = { 'At_risk': '(pali==1) | (SEMANAS<=35)',
                      'At_risk_catchup': '((pali==1) | (SEMANAS<=35))' + ' & ' + '(group == "CATCH_UP")',  
                      
                      'High_risk': '(pali==1)',
                      'High_risk_catchup': '(pali==1)' + ' & ' + '(group == "CATCH_UP")',
                      
                      'Extremely_preterms': '(SEMANAS<32 | PESO<1500)',
                      'Extremely_preterms_catchup': '(SEMANAS<32 | PESO<1500)' + ' & ' + '(group == "CATCH_UP")', #DA NA
                      
                      'Non-high-risk_preterms' : '(SEMANAS<=35) & (pali==0)',
                      'Non-high-risk_preterms_catchup' : '(SEMANAS<=35) & (pali==0)' + ' & ' + '(group == "CATCH_UP")',
                    }

filtro_comba_dict = { 'High_risk': '(pali==1) | (SEMANAS<=35)',
                    #   'High_risk_seasonal': '((pali==1) | (SEMANAS<=35))' + ' & ' + '(group == "SEASONAL")',
                      'High_risk_catchup': '((pali==1) | (SEMANAS<=35))' + ' & ' + '(group == "CATCH_UP")',  
                      
                      #'criterio_pali_normal': '(criterio_pali_normal==1)',
                      
                      'pali': '(pali==1)',
                    #   'pali_seasonal': '(pali==1)' + ' & ' + '(group == "SEASONAL")',
                      'pali_catchup': '(pali==1)' + ' & ' + '(group == "CATCH_UP")',
                      
                      'super_prematuros': '(SEMANAS<32 | PESO<1500)',
                      #'super_prematuros_seasonal': '(SEMANAS<32 | PESO<1500)' + ' & ' + '(group == "SEASONAL")', #DA NA
                      'super_prematuros_catchup': '(SEMANAS<32 | PESO<1500)' + ' & ' + '(group == "CATCH_UP")', #DA NA
                     
                      'cardiopats_1': 'RUN.isin(@lista_ruts_cardio)',
                      #'cardiopats_1_seasonal': '(RUN.isin(@lista_ruts_cardio))' + ' & ' + '(group == "SEASONAL")',
                    #   'cardiopats_1_catchup': '(RUN.isin(@lista_ruts_cardio))' + ' & ' + '(group == "CATCH_UP")', #esta wea no corre en francia
                    
                    #   'displasia': '(displa==1)',
                    #   'displasia_and_sp': '(displa==1) & (SEMANAS<32 | PESO<1500)',
                    #   'displasia_and_card': '(displa==1) & (RUN.isin(@lista_ruts_cardio))',
                    #   'displasia_and_pali': '(displa==1) & (pali==1)',
                    #   'displasia_and_pre_no_pali': '(displa==1) & (SEMANAS<=35) & (pali==0)',
                      
                      'prematuros_no_pali' : '(SEMANAS<=35) & (pali==0)',
                    #   'prematuros_no_pali_seasonal' : '(SEMANAS<=35) & (pali==0)' + ' & ' + '(group == "SEASONAL")',
                      'prematuros_no_pali_catchup' : '(SEMANAS<=35) & (pali==0)' + ' & ' + '(group == "CATCH_UP")',
                    }

def ratio_unif(ratio,dic_filtros=filtro_comba_dict):
    dic_ratio_unif = { }

    for key in dic_filtros.keys():
        dic_ratio_unif[key] = ratio
    return dic_ratio_unif

ratio_selected = '1:4'
dic_ratio_unif = ratio_unif(ratio=ratio_selected)

df_copy = df.copy()
# df_francia = df.query('(fechaIng_any.notna()) & (diag_1_leter!="J")').copy()


# list_experiments_all_born = results_match(df_case_study=df_copy,
#                                  df_control_study=df_copy, #.sample(frac=0.25, random_state=42)
#                                  filtros_dic=filtro_comba_dict,
#                                  match_vars_distance_nn=match_vars_distance_nn,
#                                  match_vars_exact_nn=match_vars_exact_nn,
#                                  match_vars_distance_IP=match_vars_distance_IP,
#                                  match_vars_exact_IP=match_vars_exact_IP_not_previ, # ['prematuro','region','group','cardio1']
#                                  weights={'edad_relativa': 1}, 
#                                  list_experiments = [],
#                                  nn=False,
#                                  ratio="1:1",
#                                  covs = ['sexo','SEMANAS','PESO'], #cardio1   ['sexo','PESO','SEMANAS'] ,'cardio1','prematuro'
#                                  dic_ratio = dic_ratio_unif) # ['sexo','SEMANAS','PESO','cardio1','prematuro']

# with open("list_experiments_all_born.pkl", "wb") as f:
#     pickle.dump(list_experiments_all_born, f)

df_copy = df.copy()
# #df_francia = df.query('(fechaIng_any.notna()) & (diag_1_leter!="J")').copy()
df_francia = df.query('(fechaIng_any.notna()) & (DIAG1.isin(@icd10_codes_fr))').copy()

list_experiments_francia_not_previ = results_match(df_case_study=df_copy,
                                 df_control_study=df_francia,
                                 filtros_dic=filtro_comba_dict,
                                 match_vars_distance_nn=match_vars_distance_nn,
                                 match_vars_exact_nn=match_vars_exact_nn,
                                 match_vars_distance_IP=match_vars_distance_IP_francia,
                                 match_vars_exact_IP= ['prematuro','Macrozones','group','PREVI'], #, #match_vars_exact_IP_not_previ, cardio1
                                 weights=weights,
                                 list_experiments = [],
                                 nn=False,
                                 ratio="1:1",
                                 covs =  ['sexo','SEMANAS','PESO'])  #cardio1 



with open("list_experiments_francia_not_previ.pkl", "wb") as f:
    pickle.dump(list_experiments_francia_not_previ, f)

df_final = summary_eff_aux(list_experiments_francia_not_previ = list_experiments_francia_not_previ)
# df_final = summary_eff_aux(list_experiments_all_born=list_experiments_all_born) #list_experiments_all_born, list_experiments_francia_not_previ
#print(df_final)

with pd.ExcelWriter(path_data/"resultados_comparados_copy_nuevo.xlsx", engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    df_final.to_excel(writer, sheet_name=f"T2_francia_1a{ratio_selected.split(':')[1]}_rev_final") #Tabla_2_ratio_1_1_unif_upc   Tabla_2_1a7_francia_oficial
    

# df_problemas = df_problemas_matching(list_1=list_experiments_all_born,
                                    # list_2=list_experiments_francia_not_previ)
# with pd.ExcelWriter(path_data/"resultados_comparados_copy.xlsx", engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
#     df_problemas.to_excel(writer, sheet_name="df_problemas_ratio_optim")


############################################################# UPC ##############################################################

# filtro_comba_dict = {'UPC_general': 'MARCA==0'}

# df_copy = df_upc.copy()

# list_experiments_all_born = results_match(df_case_study=df_copy,
#                                  df_control_study=df_copy.sample(frac=0.10, random_state=42),
#                                  filtros_dic=filtro_comba_dict,
#                                  match_vars_distance_nn=match_vars_distance_nn,
#                                  match_vars_exact_nn=match_vars_exact_nn,
#                                  match_vars_distance_IP=match_vars_distance_IP,
#                                  match_vars_exact_IP=match_vars_exact_IP_not_previ,
#                                  weights=weights, 
#                                  list_experiments = list_experiments_all_born, #list_experiments_all_born
#                                  nn=False,
#                                  ratio="1:1",
#                                  covs = ['sexo','PESO','SEMANAS'], #displa, cardio1
#                                  dic_ratio = dic_ratio_unif) # ['sexo','SEMANAS','PESO','cardio1','prematuro']


# df_copy = df_upc.copy()
# # df_francia = df.query('(fechaIng_any.notna()) & (diag_1_leter!="J")').copy()
# df_francia = df_upc.query('(fechaIng_any.notna()) & (DIAG1.isin(@icd10_codes_fr))').copy()

# list_experiments_francia_not_previ = results_match(df_case_study=df_copy,
#                                  df_control_study=df_francia,
#                                  filtros_dic=filtro_comba_dict,
#                                  match_vars_distance_nn=match_vars_distance_nn,
#                                  match_vars_exact_nn=match_vars_exact_nn,
#                                  match_vars_distance_IP=match_vars_distance_IP_francia,
#                                  match_vars_exact_IP=match_vars_exact_IP_not_previ,
#                                  weights=weights, 
#                                  list_experiments = list_experiments_francia_not_previ, #list_experiments_francia_not_previ
#                                  nn=False,
#                                  ratio="1:1",
#                                  covs =  ['sexo','SEMANAS','PESO','PREVI']) #cardio1

# # df_copy = df.copy()
# #list_experiments_all_born_upc = list_experiments_all_born.copy()
# list_experiments_francia_not_previ_upc = list_experiments_francia_not_previ.copy()

# # with open("list_experiments_all_born.pkl", "wb") as f:
# #     pickle.dump(list_experiments_all_born_upc, f)
# with open("list_experiments_francia_not_previ.pkl", "wb") as f:
#     pickle.dump(list_experiments_francia_not_previ_upc, f)

# df_final = summary_eff_aux(list_experiments_all_born, list_experiments_francia_not_previ) #list_experiments_all_born,
# print(df_final)

# # df_problemas = df_problemas_matching(#list_1=list_experiments_all_born,
# #                                     list_2=list_experiments_francia_not_previ)


# ################################################## SAVEING ##########################################
    
# with pd.ExcelWriter(path_data/"resultados_comparados_copy.xlsx", engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
#     df_final.to_excel(writer, sheet_name="Tabla_palis_upc") #Tabla_2_ratio_1_1_unif_upc