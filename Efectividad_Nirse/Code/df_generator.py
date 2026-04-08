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
from lifelines import KaplanMeierFitter
import warnings
from preproces_prod3 import *
from scipy.stats import norm
import gspread
from oauth2client.service_account import ServiceAccountCredentials
warnings.filterwarnings("ignore")

path_actual = Path.cwd()
print(path_actual)

path_data = path_actual.parent / 'Data'
print(path_data)

print('start')
df_cox_vrs_Dall, df_cox_upc_Dall, df_f_vrs_Dall, df_f_upc_Dall = call_data_cox('NAC_RNI_EGRESOS_ENTREGA_ISCI_20_12_2024_encr.csv',40,group_age=True,weeks_inm=False)

print('here1')
df_cox_lrti_Dall, df_f_LRTI_Dall = call_data_cox_auxiliar('NAC_RNI_EGRESOS_ENTREGA_ISCI_20_12_2024_encr.csv',40,group_age=True,fecha_vrs='fechaIng_LRTI',lrti_name='LRTI_Flag')

print('here2')
df_cox_any_Dall, df_f_any_Dall = call_data_cox_auxiliar('NAC_RNI_EGRESOS_ENTREGA_ISCI_20_12_2024_encr.csv',40,group_age=True,fecha_vrs='fechaIng_any')

print('here3')
df_cox_vrs_Dall['group_age'] = df_cox_vrs_Dall[['si_1_meses', 'si_2_meses', 'si_3_meses', 'si_4_meses', 'si_5_meses', 'si_6_meses']].apply(
    lambda row: row.index[row == 1].max(), axis=1).fillna('si_0_meses')

print('here4')
df_cox_upc_Dall['group_age'] = df_cox_upc_Dall[['si_1_meses', 'si_2_meses', 'si_3_meses', 'si_4_meses', 'si_5_meses', 'si_6_meses']].apply(
    lambda row: row.index[row == 1].max(), axis=1).fillna('si_0_meses')

print('here5')
df_cox_lrti_Dall['group_age'] = df_cox_lrti_Dall[['si_1_meses', 'si_2_meses', 'si_3_meses', 'si_4_meses', 'si_5_meses', 'si_6_meses']].apply(
    lambda row: row.index[row == 1].max(), axis=1).fillna('si_0_meses')

print('here6')
df_cox_any_Dall['group_age'] = df_cox_any_Dall[['si_1_meses', 'si_2_meses', 'si_3_meses', 'si_4_meses', 'si_5_meses', 'si_6_meses']].apply(
    lambda row: row.index[row == 1].max(), axis=1).fillna('si_0_meses')


print('here7')
df_cox_vrs_Dall.to_csv(path_data/'df_vrs_s39_2012_all_meses_Dall.csv', index= False)

print('here8')
df_cox_upc_Dall.to_csv(path_data/'df_upc_s39_2012_all_meses_Dall.csv', index= False)

print('here9')
df_cox_lrti_Dall.to_csv(path_data/'df_lrti_s39_2012_all_meses_Dall.csv', index= False)

print('here10')
df_cox_any_Dall.to_csv(path_data/'df_any_s39_2012_all_meses_Dall.csv', index= False)

print('end')

