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
path_data = path_actual.parent / 'Data'

df_pf = pre_filtred(df_name='NAC_RNI_EGRESOS_ENTREGA_ISCI_20_12_2024_encr.csv',lrti_name='LRTI_Flag')
df_f_any, df_f_LRTI, df_f_vrs, df_f_upc = filtros_IH_new(df_pf)

print('que')
df_vrs = pd.read_csv(path_data/'df_vrs_s39_2012_all_meses_Dall.csv')
df_upc = pd.read_csv(path_data/'df_upc_s39_2012_all_meses_Dall.csv')
df_LRTI = pd.read_csv(path_data/'df_LRTI_s39_2012_all_meses_Dall.csv')
df_any = pd.read_csv(path_data/'df_any_s39_2012_all_meses_Dall.csv')

print('que2')
covs_base = ['start', 'inmunizado', 'stop', 'RUN', 'sexo','region', 'group_age','SEMANAS']

events_and_dataframes = {
    'event_vrs': df_vrs,
    'event_any': df_any,
    'event_LRTI': df_LRTI,
    'event_upc': df_upc
}

cox_results = {}

print('que3')
for event, df_cox in events_and_dataframes.items():
    covs = covs_base + [event]  # Crear la lista de covariables dinámicamente
    cox_results[event] = cox_return(df_cox, covs=covs, prematuros=False).reset_index()
print('que4')
cox_vrs = cox_results['event_vrs']
cox_any = cox_results['event_any']
cox_LRTI = cox_results['event_LRTI']
cox_upc = cox_results['event_upc']

T_final_tabla = pd.to_datetime('2024-09-30')

df_f_any['anyIng_time_days'] = (T_final_tabla - df_f_any[['fechaIng_any', 'fecha_nac']].max(axis=1, skipna=True)).dt.days
df_f_LRTI['LRTI_time_days'] = (T_final_tabla - df_f_LRTI[['fechaIng_LRTI', 'fecha_nac']].max(axis=1, skipna=True)).dt.days
df_f_vrs['VRS_time_days'] = (T_final_tabla - df_f_vrs[['fechaIng_vrs', 'fecha_nac']].max(axis=1, skipna=True)).dt.days
df_f_upc['upc_time_days'] = (T_final_tabla - df_f_upc[['fecha_upc_vrs', 'fecha_nac']].max(axis=1, skipna=True)).dt.days

df_f = df_f_vrs.copy()
coef_inmunizado = cox_vrs.query('covariate=="inmunizado"').effectiveness.values[0]

lower_bound = cox_vrs.query('covariate=="inmunizado"').eff_lower_95.values[0]
upper_bound =  cox_vrs.query('covariate=="inmunizado"').eff_upper_95.values[0]

Effective = f"{coef_inmunizado:.2f} ({lower_bound:.2f}-{upper_bound:.2f})"

VRS_event_nonrec = df_f[df_f['inmunizado'] == 0]['event_vrs'].sum()
VRS_event_rec = df_f[df_f['inmunizado'] == 1]['event_vrs'].sum()

VRS_person_year_nonrec = df_f[df_f['inmunizado'] == 0]['VRS_time_days'].sum()
VRS_person_year_rec = df_f[df_f['inmunizado'] == 1]['VRS_time_days'].sum()

df_f = df_f_upc.copy()
coef_inmunizado_upc = cox_upc.query('covariate=="inmunizado"').effectiveness.values[0]

lower_bound_upc = cox_upc.query('covariate=="inmunizado"').eff_lower_95.values[0]
upper_bound_upc = cox_upc.query('covariate=="inmunizado"').eff_upper_95.values[0]

Effective_upc = f"{coef_inmunizado_upc:.2f} ({lower_bound_upc:.2f}-{upper_bound_upc:.2f})"

VRS_event_nonrec_upc = df_f[df_f['inmunizado'] == 0]['event_upc'].sum()
VRS_event_rec_upc = df_f[df_f['inmunizado'] == 1]['event_upc'].sum()

VRS_person_year_nonrec_upc = df_f[df_f['inmunizado'] == 0]['upc_time_days'].sum()
VRS_person_year_rec_upc = df_f[df_f['inmunizado'] == 1]['upc_time_days'].sum()

df_f = df_f_LRTI.copy()
coef_inmunizado_LRTI = cox_LRTI.query('covariate=="inmunizado"').effectiveness.values[0]

lower_bound_LRTI = cox_LRTI.query('covariate=="inmunizado"').eff_lower_95.values[0]
upper_bound_LRTI = cox_LRTI.query('covariate=="inmunizado"').eff_upper_95.values[0]

Effective_LRTI = f"{coef_inmunizado_LRTI:.2f} ({lower_bound_LRTI:.2f}-{upper_bound_LRTI:.2f})"

LRTI_event_nonrec = df_f[df_f['inmunizado'] == 0]['event_LRTI'].sum()
LRTI_event_rec = df_f[df_f['inmunizado'] == 1]['event_LRTI'].sum()

LRTI_person_year_nonrec = df_f[df_f['inmunizado'] == 0]['LRTI_time_days'].sum()
LRTI_person_year_rec = df_f[df_f['inmunizado'] == 1]['LRTI_time_days'].sum()

df_f = df_f_any.copy()
coef_inmunizado_any = cox_any.query('covariate=="inmunizado"').effectiveness.values[0]

lower_bound_any = cox_any.query('covariate=="inmunizado"').eff_lower_95.values[0]
upper_bound_any = cox_any.query('covariate=="inmunizado"').eff_upper_95.values[0]

Effective_any = f"{coef_inmunizado_any:.2f} ({lower_bound_any:.2f}-{upper_bound_any:.2f})"

any_event_nonrec = df_f[df_f['inmunizado'] == 0]['event_any'].sum()
any_event_rec = df_f[df_f['inmunizado'] == 1]['event_any'].sum()

any_person_year_nonrec = df_f[df_f['inmunizado'] == 0]['anyIng_time_days'].sum()
any_person_year_rec = df_f[df_f['inmunizado'] == 1]['anyIng_time_days'].sum()

df_f = df_f_vrs.copy()
total_overall = len(df_f)
total_non_recipients = len(df_f[df_f['inmunizado'] == 0])
total_recipients = len(df_f[df_f['inmunizado'] == 1])


nirse_index = f"Nirse recipients \n(N={total_recipients})"
non_nirse_index = f"Nirse non recipients \n(N={total_non_recipients})"

column_tuples2 = [
    (non_nirse_index, "Events"),
    ("", "Person-years"),
    (nirse_index, "Events"),
    (" ", "Person-years"),
    ("Hazard_rate (95% CI)", ""),
    ("Effectiveness (95% CI)", "")
]

columns2 = pd.MultiIndex.from_tuples(column_tuples2)

row_index2 = pd.Index([
    "RSV-related LRTI hospitalisation",
    "Severe RSV-related LRTI with intensive care unit admission",
    "All-cause LRTI hospitalisation",
    "All-cause hospitalisation"
])


df_table2 = pd.DataFrame(index=row_index2, columns=columns2)

df_table2.loc["RSV-related LRTI hospitalisation", (non_nirse_index, "Events")] = VRS_event_nonrec
df_table2.loc["RSV-related LRTI hospitalisation", ("", "Person-years")] = round(VRS_person_year_nonrec/365,2)
df_table2.loc["RSV-related LRTI hospitalisation", (nirse_index, "Events")] = VRS_event_rec
df_table2.loc["RSV-related LRTI hospitalisation", (" ", "Person-years")]= round(VRS_person_year_rec/365,2)
df_table2.loc["RSV-related LRTI hospitalisation", ("Hazard_rate (95% CI)", "")] = pd.NA
df_table2.loc["RSV-related LRTI hospitalisation", ("Effectiveness (95% CI)", "")] = Effective

df_table2.loc["Severe RSV-related LRTI with intensive care unit admission", (non_nirse_index, "Events")] = VRS_event_nonrec_upc
df_table2.loc["Severe RSV-related LRTI with intensive care unit admission", ("", "Person-years")] = round(VRS_person_year_nonrec_upc/365,2)
df_table2.loc["Severe RSV-related LRTI with intensive care unit admission", (nirse_index, "Events")] = VRS_event_rec_upc
df_table2.loc["Severe RSV-related LRTI with intensive care unit admission", (" ", "Person-years")]= round(VRS_person_year_rec_upc/365,2)
df_table2.loc["Severe RSV-related LRTI with intensive care unit admission", ("Hazard_rate (95% CI)", "")] = pd.NA
df_table2.loc["Severe RSV-related LRTI with intensive care unit admission", ("Effectiveness (95% CI)", "")] = Effective_upc

df_table2.loc["All-cause LRTI hospitalisation", (non_nirse_index, "Events")] = LRTI_event_nonrec
df_table2.loc["All-cause LRTI hospitalisation", ("", "Person-years")] = round(LRTI_person_year_nonrec/365,2)
df_table2.loc["All-cause LRTI hospitalisation", (nirse_index, "Events")] = LRTI_event_rec
df_table2.loc["All-cause LRTI hospitalisation", (" ", "Person-years")]= round(LRTI_person_year_rec/365,2)
df_table2.loc["All-cause LRTI hospitalisation", ("Hazard_rate (95% CI)", "")] = pd.NA
df_table2.loc["All-cause LRTI hospitalisation", ("Effectiveness (95% CI)", "")] = Effective_LRTI

df_table2.loc["All-cause hospitalisation", (non_nirse_index, "Events")] = any_event_nonrec
df_table2.loc["All-cause hospitalisation", ("", "Person-years")] = round(any_person_year_nonrec/365,2)
df_table2.loc["All-cause hospitalisation", (nirse_index, "Events")] = any_event_rec
df_table2.loc["All-cause hospitalisation", (" ", "Person-years")]= round(any_person_year_rec/365,2)
df_table2.loc["All-cause hospitalisation", ("Hazard_rate (95% CI)", "")] = pd.NA
df_table2.loc["All-cause hospitalisation", ("Effectiveness (95% CI)", "")] = Effective_any


display(df_table2)