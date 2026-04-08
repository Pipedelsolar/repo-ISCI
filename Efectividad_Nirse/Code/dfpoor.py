import pandas as pd
import math
from datetime import datetime, timedelta
from difflib import SequenceMatcher
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import os
from pathlib import Path
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import seaborn as sns
import numpy as np
from lifelines.utils import to_long_format
from lifelines.utils import add_covariate_to_timeline
from lifelines import CoxTimeVaryingFitter
from lifelines import CoxPHFitter
import warnings
from Preproces_prod2 import *
from patsy import dmatrices
warnings.filterwarnings("ignore")


path_actual = Path.cwd()
path_data = path_actual.parent / 'Data'

df_indic = pd.read_excel(path_data/'Estimaciones_Indice_Pobreza_Multidimensional_Comunas_2022.xlsx', skiprows=2) 
df_tas = pd.read_excel(path_data/'Estimaciones_Tasa_Pobreza_Ingresos_Comunas_2022.xlsx', skiprows=2) 

df_tasa = df_tas.head(345)
df_tasa = df_tasa[['Código',
                   'Porcentaje de personas en situación de pobreza por ingresos 2022']].rename(columns={'Código': 'COMUNA', 'Porcentaje de personas en situación de pobreza por ingresos 2022': 'percent_poor'})

df_indice = df_indic.head(345)
df_indice = df_indice[['Código','Porcentaje de personas en situación de pobreza multidimensional 2022']].rename(columns={'Código': 'COMUNA', 'Porcentaje de personas en situación de pobreza multidimensional 2022': 'percent_poor_multidim'})

df_com_poor = df_indice.merge(df_tasa,on='COMUNA', how='left')
df_com_poor.to_csv(path_data/'df_com_poor.csv', index=False)