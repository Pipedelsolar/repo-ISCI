import pandas as pd
from datetime import datetime, timedelta
from difflib import SequenceMatcher
import matplotlib.pyplot as plt
import pickle
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
path_nirsecl = path_actual.parent.parent/'Nirse_cl' / 'Data'

with open(path_data/'lista_ruts_cardio.pkl', 'rb') as f:
    lista_ruts_cardio = pickle.load(f)

with open(path_data/'lista_ruts_preterms.pkl', 'rb') as f:
    lista_ruts_preterms = pickle.load(f)
    
df_pf = pre_filtred(df_name='NAC_RNI_EGRESOS_ENTREGA_ISCI_11_04_2025_encr.csv',lrti_name='LRTI_Flag')
#_, _, df_f_vrs, _ = filtros_IH_new(df_pf)

path_nirsecl = path_actual.parent.parent/'Nirse_cl' / 'Data'
df_historic = pd.read_csv(path_nirsecl/"data_sintrib.csv").drop(columns=['Unnamed: 0'])  #data.csv #.query('RUN.isin(@ruts_cariopaticos_1.RUN.unique())')

filtro_cardio = 'RUN.isin(@lista_ruts_cardio)'
filtro_prematuro = 'RUN.isin(@lista_ruts_preterms)'
filtro_any = '(RUN.isin(@lista_ruts_cardio)) | (RUN.isin(@lista_ruts_preterms))'

df_historic_cardio = df_historic.query(filtro_cardio).reset_index(drop=True)
df_historic_prematuro = df_historic.query(filtro_prematuro).reset_index(drop=True)

df_nac_hist = pd.read_csv(path_data / "NAC_DEF_EGRE3.csv", encoding='latin1', sep=';')
pf_nac = pre_filtred_nac_hist(df_name='NAC_DEF_EGRE3.csv')

pf_all = pd.concat([pf_nac,df_pf.assign(FECHA_NAC = lambda x: x.fecha_nac)]).reset_index(drop=True)
_, _, nac_filtered_vrs, _ = filtros_IH_nac(pf_all)

df_pf.to_csv(path_data / "df_pf.csv", index=False)
df_historic.to_csv(path_data / "df_historic.csv", index=False)
df_historic_cardio.to_csv(path_data / "df_historic_cardio.csv", index=False)
df_historic_prematuro.to_csv(path_data / "df_historic_prematuro.csv", index=False)

pf_nac.to_csv(path_data / "pf_nac.csv", index=False)
pf_all.to_csv(path_data / "pf_all.csv", index=False)
nac_filtered_vrs.to_csv(path_data / "nac_filtered_vrs.csv", index=False)

df_nac_per_year_filt_cardio = count_nacs_filter(pf_all, nac_filtered_vrs, filtro_cardio)
df_nac_per_year_filt_prematuro = count_nacs_filter(pf_all, nac_filtered_vrs, filtro_prematuro)
df_nac_per_year_filt_any = count_nacs_filter(pf_all, nac_filtered_vrs, filtro_any)

nacimientos_per_year_cardio = df_nac_per_year_filt_cardio[['elegibilidad_alt','nacs_filtred']].rename(columns={'nacs_filtred':'nacimientos','elegibilidad_alt':'year'})
nacimientos_per_year_prematuro = df_nac_per_year_filt_prematuro[['elegibilidad_alt','nacs_filtred']].rename(columns={'nacs_filtred':'nacimientos','elegibilidad_alt':'year'})
nacimientos_per_year_any = df_nac_per_year_filt_any[['elegibilidad_alt','nacs_filtred']].rename(columns={'nacs_filtred':'nacimientos','elegibilidad_alt':'year'})

nacimientos_per_year_cardio.to_csv(path_data / "nacimientos_per_year_cardio.csv", index=False)
nacimientos_per_year_prematuro.to_csv(path_data / "nacimientos_per_year_prematuro.csv", index=False)
nacimientos_per_year_any.to_csv(path_data / "nacimientos_per_year_any.csv", index=False)

df_pf.to_csv(path_nirsecl / "df_pf.csv", index=False)
df_historic.to_csv(path_nirsecl / "df_historic.csv", index=False)
df_historic_cardio.to_csv(path_nirsecl / "df_historic_cardio.csv", index=False)
df_historic_prematuro.to_csv(path_nirsecl / "df_historic_prematuro.csv", index=False)
pf_nac.to_csv(path_nirsecl / "pf_nac.csv", index=False)
pf_all.to_csv(path_nirsecl / "pf_all.csv", index=False)
nac_filtered_vrs.to_csv(path_nirsecl / "nac_filtered_vrs.csv", index=False)

nacimientos_per_year_cardio.to_csv(path_nirsecl / "nacimientos_per_year_cardio.csv", index=False)
nacimientos_per_year_prematuro.to_csv(path_nirsecl / "nacimientos_per_year_prematuro.csv", index=False)
nacimientos_per_year_any.to_csv(path_nirsecl / "nacimientos_per_year_any.csv", index=False)



