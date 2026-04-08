import pandas as pd
from datetime import datetime, timedelta
from difflib import SequenceMatcher
import matplotlib.pyplot as plt
import os
import numpy as np
from lifelines.utils import to_long_format
from lifelines.utils import add_covariate_to_timeline
from scipy.stats import ks_2samp
import warnings
from preproces_prod3 import *
from matching_case_control import *

warnings.filterwarnings("ignore")

import pickle

with open(path_data/'lista_ruts_cardio.pkl', 'rb') as f:
    lista_ruts_cardio = pickle.load(f)

#df_vrs_match_case, _ = call_data('NAC_RNI_EGRESOS_ENTREGA_ISCI_11_04_2025_encr.csv')
df_pf = (
        pre_filtred(df_name='NAC_RNI_EGRESOS_ENTREGA_ISCI_11_04_2025_encr.csv')
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
          #  Macrozones = lambda x: x.region.map(region_to_macrozone_agencia).replace(leo_zonas_rename),
            super_preterm = lambda x: ((x.SEMANAS<32) | (x.PESO<1500)).astype(int),
        )
        .copy()
    )

df_pf_sd = (
    df_pf
    .copy()
    .assign(fechaIng_any_modify = lambda x: x.fechaIng_any,
                        fechaIng_vrs_modify = lambda x: x.fechaIng_vrs,
                        fecha_def = lambda x: pd.to_datetime(x['FECHA_DEFUNCION'], format='%d%b%Y'),
                        muerto_pre_campain = lambda x: x.fecha_def < pd.to_datetime("2024-04-01"),
                        muerto_7days_birth = lambda x: x.fecha_def <= x.fecha_nac + pd.DateOffset(days=7))
    # .query('muerto_pre_campain==False') 
    # .query('muerto_7days_birth==False')
    #.sort_values(by='fechaIng_vrs').drop_duplicates(subset=['RUN'], keep='first')
    #.query('VIVO=="SI"')
    .query('((SEMANAS<32) | (PESO<1500) | (RUN.isin(@lista_ruts_cardio))) | ((32<=SEMANAS<=34) & (PESO>2500)) | ((PESO>=1500) & (SEMANAS==35)) & ~(RUN.isin(@lista_ruts_cardio))')
    )


df_pf_sd.loc[df_pf_sd.fechaIng_any_modify <= df_pf_sd.fecha_nac + pd.DateOffset(days=7), ['fechaIng_any_modify','fechaIng_vrs_modify','fecha_upc_vrs','fechaIng_LRTI']] = pd.NaT

filters = {
    # 'G.a.w > 42': (df_pf_sd.SEMANAS > 42),
    # 'G.a.w < 24': (df_pf_sd.SEMANAS < 24),
    # 'G.weight <=  p_l': (df_pf_sd.PESO <= df_pf_sd.p_00001_lognormal),
    # 'G.weight >= p_u': (df_pf_sd.PESO >= df_pf_sd.p_99999_lognormal),
    'death_at_born': (df_pf_sd['VIVO']=="NO") & (df_pf_sd['muerto_7days_birth']),
    'curve_peso_edadGest': ~((df_pf_sd['SEMANAS']>=24) & (df_pf_sd['SEMANAS']<=42) & (df_pf_sd['PESO'] > 500) & (df_pf_sd['PESO'] <= df_pf_sd['p_99999_lognormal'])),#((df_pf_sd.PESO < df_pf_sd.p_00001_lognormal) | (df_pf_sd.PESO > df_pf_sd.p_99999_lognormal) | (df_pf_sd.SEMANAS > 42) | (df_pf_sd.SEMANAS < 24)),
    'Intersex': (df_pf_sd.SEXO == 9),
   # 'Unidentifiable ID': (df_pf_sd.MARCA==1),
   # 'Mom age': ((df_pf_sd.EDAD_M > 51) | (df_pf_sd.EDAD_M < 12)),
   # 'vrs_7_days_birth': (df_pf_sd.RUN.isin(df_pf_sd[df_pf_sd.fechaIng_any >= df_pf_sd.fecha_nac].RUN.unique())) & (df_pf_sd.RUN.isin(df_pf_sd[(df_pf_sd.fechaIng_vrs <= df_pf_sd.fecha_nac + pd.DateOffset(days=7))].RUN.unique())), #(df_pf_sd.fechaIng_vrs <= df_pf_sd.fecha_nac + pd.DateOffset(days=7)), #
    'Inconsistent dates': (df_pf_sd.RUN.isin(df_pf_sd[df_pf_sd.fechaIng_any < df_pf_sd.fecha_nac].RUN.unique())), #
    'VRS Hospitalization date before start of campaign': (df_pf_sd.RUN.isin(df_pf_sd[df_pf_sd.fechaIng_vrs_modify < pd.to_datetime("2024-04-01")].RUN.unique())),
}

# Crear la tabla vacía con los nombres de los filtros
filter_names = list(filters.keys())
intersection_table = pd.DataFrame(index=filter_names, columns=filter_names + ['Total'])

df_diagonal = {}

# Rellenar la tabla con las intersecciones
for f1 in filter_names:
    for f2 in filter_names:
        if f1 == f2:
            super_filter = filters[f1]
            for f3 in filter_names:
                if f3 == f1:
                    continue
                else:
                    super_filter = super_filter & (~filters[f3])
                    
                
            count = df_pf_sd[super_filter].RUN.nunique() #len(df_pf_sd[super_filter])
            df_diagonal[f1] = count
        else:
            # Intersección: cuántos datos elimina al aplicar f1 y f2 simultáneamente
            combined_filter = filters[f1] & filters[f2]
            count = df_pf_sd[combined_filter].RUN.nunique() #len(df_pf_sd[combined_filter])
        
        intersection_table.loc[f1, f2] = count

    # Calcular el total general (datos eliminados al aplicar f1 sin exclusividad)
    intersection_table.loc[f1, 'Total'] = df_pf_sd[filters[f1]].RUN.nunique() - intersection_table.loc[f1, f1]  #len(df_pf_sd[filters[f1]])


# diag_serie = pd.Series(df_diagonal)
# total_compare = intersection_table.Total - diag_serie

with pd.ExcelWriter(path_data/"Tabla_exclusion.xlsx", engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    intersection_table.to_excel(writer, sheet_name="Tabla_exclusion")



# df_pf_sd = df_pf.assign(fechaIng_any_modify = lambda x: x.fechaIng_any,
#                         fechaIng_vrs_modify = lambda x: x.fechaIng_vrs)

# 1) Crear columnas booleanas en el DataFrame para cada filtro
for fname, cond in filters.items():
    df_pf_sd[fname] = cond  # Esto creará o sobrescribirá la columna con valores True/False

filter_names = list(filters.keys())
# df_pf_sd['num_filters_cumplidos'] = df_pf_sd[filter_names].sum(axis=1)

ruts_filts = df_pf_sd.assign(num_filters_cumplidos = lambda x: x[filter_names].sum(axis=1)).groupby('RUN',as_index=False).agg({'num_filters_cumplidos': 'max'})

df_pf_sd = df_pf_sd.merge(ruts_filts, on='RUN', how='left')

# 3) Construir la tabla de intersecciones para “3 o más” filtros
intersection_table_3plus = pd.DataFrame(
    index=filter_names,
    columns=filter_names + ['Total']
)

for f1 in filter_names:

    intersection_table_3plus.loc[f1, 'Total'] = df_pf_sd[df_pf_sd[f1] & (df_pf_sd['num_filters_cumplidos'] >= 3)].RUN.nunique()
    
    for f2 in filter_names:
        if f1 == f2:
            # Podrías dejar vacío, o poner 0, o un guión
            intersection_table_3plus.loc[f1, f2] = 'x'
        else:
            # Marcamos “x” si existe AL MENOS un registro que cumple:
            #   f1 == True, f2 == True y num_filters_cumplidos >= 3
            mask_f1_f2_3plus = (
                df_pf_sd[f1] & 
                df_pf_sd[f2] & 
                (df_pf_sd['num_filters_cumplidos'] >= 3)
            )
            intersection_table_3plus.loc[f1, f2] = 'x' if df_pf_sd[mask_f1_f2_3plus].RUN.nunique()>=1 else ''
            
df_diagonal_real = pd.DataFrame([df_diagonal])
df_diagonal_real.index = ['Diagonal']
intersection_table_3plus = pd.concat([intersection_table_3plus, df_diagonal_real], ignore_index=False)

with pd.ExcelWriter(path_data/"Tabla_exclusion.xlsx", engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    intersection_table_3plus.to_excel(writer, sheet_name="Tabla_exclusion_intersect_3")
    