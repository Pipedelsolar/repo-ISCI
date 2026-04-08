
import os 
import sys

import pandas as pd
import numpy as np

from preproces_prod3 import *
    
local_path = os.getcwd()
code_root = os.path.abspath(os.path.join(local_path, '..', 'Code'))

if code_root not in sys.path:
    sys.path.insert(0, code_root)
code_root = os.path.abspath(os.path.join(local_path, '../..', 'inv'))

if code_root not in sys.path:
    sys.path.insert(0, code_root)

from matching_case_control import call_data,analyze_vrs_data,match_nn_max_dist_weigths,comparar_medias_test, charly_mip, charly_double_mip, mylogit, results_match, tabla_marcel, tabla_final,summary_eff,tabla_1_match

import inv

from IPython.core.display import display
display.max_output_lines = 500  # Adjust the number as needed
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

np.random.seed(42)

path_actual = Path.cwd()
path_data = path_actual.parent / 'Data'

# df_vrs_match_case = pd.read_csv(path_data/'df_vrs_match_case_11_04_25.csv', index_col=0) 
# df_upc_match_case = pd.read_csv(path_data/'df_upc_match_case_11_04_25.csv', index_col=0)

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


def _tabla_1_single(
    df,
    id_col='RUN',
    group_col='Group',
    case_label='Caso',
    control_label='Control',
    decimales=2
):
    import numpy as np
    import pandas as pd

    df = df.copy()
    df['grupo_tabla'] = np.where(df[id_col] == df[group_col], case_label, control_label)
    df['sexo'] = df['sexo'].map({0: 'Female', 1: 'Male'})

    variables_numericas_mediana = ['edad_relativa']
    variables_numericas_media = ['SEMANAS', 'PESO']
    variables_categoricas = ['sexo']
    #macrozona_var = 'Macrozones'

    etiquetas_variables = {
        'edad_relativa': 'Age (days)',
        'SEMANAS': 'Gestational age at birth (weeks)',
        'PESO': 'Birth weight (grams)',
        'sexo': 'Sex',
    #    'Macrozones': 'Macro-zone',
        'cardio1': 'Congenital heart disease',
        'prematuro': 'Prematurity',
        'event_upc': 'PICU admission',
        'super_preterm':'Prematurity Extreme'
    }

    grupos = df['grupo_tabla'].unique()
    tabla = []
    index = []
    n_por_grupo = {g: len(df[df['grupo_tabla'] == g]) for g in grupos}
    nombres_columnas = {g: f"{g} (n = {n_por_grupo[g]})" for g in grupos}

    # ---- Variables numéricas con mediana (IQR)
    for var in variables_numericas_mediana:
        fila = {}
        for g in grupos:
            subset = df[df['grupo_tabla'] == g][var]
            mediana = subset.median()
            iqr = subset.quantile(0.75).round(1).astype(str) + '-' + subset.quantile(0.25).round(1).astype(str)
            # fila[nombres_columnas[g]] = f'{mediana:.{decimales}f} ({iqr:.{decimales}f})'
            fila[nombres_columnas[g]] = f'{mediana:.{decimales}f} (' + iqr + ')'
        index.append((etiquetas_variables.get(var, var), ''))
        tabla.append(fila)

        # Rango por edad relativa: <90, 90-180, >180
        bins = [0, 90, 180, float('inf')]
        labels = ['<90', '90-180', '>180']
        df['edad_intervalo'] = pd.cut(df['edad_relativa'], bins=bins, labels=labels, right=False)
        for val in labels:
            fila = {}
            for g in grupos:
                subset = df[df['grupo_tabla'] == g]
                total = len(subset)
                cuenta = (subset['edad_intervalo'] == val).sum()
                porcentaje = 100 * cuenta / total if total > 0 else 0
                fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
            index.append((etiquetas_variables.get(var, var), val))
            tabla.append(fila)

    # ---- Variables categóricas comunes
    for var in variables_categoricas:
        index.append((etiquetas_variables.get(var, var), ''))
        tabla.append({})  # Fila vacía para separar
        niveles = df[var].dropna().unique()
        for val in niveles:
            fila = {}
            for g in grupos:
                subset = df[df['grupo_tabla'] == g]
                total = len(subset)
                cuenta = (subset[var] == val).sum()
                porcentaje = 100 * cuenta / total if total > 0 else 0
                fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
            index.append((etiquetas_variables.get(var, var), val))
            tabla.append(fila)

    # ---- Variables numéricas con media (SD)
    for var in variables_numericas_media:
        fila = {}
        for g in grupos:
            subset = df[df['grupo_tabla'] == g][var]
            media = subset.mean()
            sd = subset.std()
            fila[nombres_columnas[g]] = f'{media:.{decimales}f} ({sd:.{decimales}f})'
        index.append((etiquetas_variables.get(var, var), ''))
        tabla.append(fila)

        if var == 'SEMANAS':
            # Rango para SEMANAS: 24-28 y 29-35
            bins = [0, 29, 36]  # hasta <29 y <36
            labels = ['24-28', '29-35']
            df['semana_intervalo'] = pd.cut(df['SEMANAS'], bins=bins, labels=labels, right=False)
            for val in labels:
                fila = {}
                for g in grupos:
                    subset = df[df['grupo_tabla'] == g]
                    total = len(subset)
                    cuenta = (subset['semana_intervalo'] == val).sum()
                    porcentaje = 100 * cuenta / total if total > 0 else 0
                    fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
                index.append((etiquetas_variables.get(var, var), val))
                tabla.append(fila)
    index.append(('Risk groups', ''))
    tabla.append({})
    for var in ['cardio1', 'prematuro','super_preterm']:
        fila = {}
        for g in grupos:
            subset = df[df['grupo_tabla'] == g]
            total = len(subset)
            cuenta = (subset[var] == 1).sum()
            porcentaje = 100 * cuenta / total if total > 0 else 0
            fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
        index.append(('Risk groups', etiquetas_variables.get(var, var)))
        tabla.append(fila)

    # ---- Macrozonas
    # index.append((etiquetas_variables.get(macrozona_var, macrozona_var), ''))
    # tabla.append({})
    # niveles = df[macrozona_var].dropna().unique()
    # for val in niveles:
    #     fila = {}
    #     for g in grupos:
    #         subset = df[df['grupo_tabla'] == g]
    #         total = len(subset)
    #         cuenta = (subset[macrozona_var] == val).sum()
    #         porcentaje = 100 * cuenta / total if total > 0 else 0
    #         fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
    #     index.append((etiquetas_variables.get(macrozona_var, macrozona_var), val))
    #     tabla.append(fila)

    # ---- PICU
    fila = {}
    for g in grupos:
        subset = df[df['grupo_tabla'] == g]
        total = len(subset)
        cuenta = (subset['cama'] == 'UPC').sum()
        porcentaje = 100 * cuenta / total if total > 0 else 0
        fila[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
    index.append((etiquetas_variables.get('cama', 'PICU'), ''))
    tabla.append(fila)

    # -------------------------------------------------------------------------
    # A) Nirsevimab Recipients (inmunizado=1 => "Yes", inmunizado=0 => "No")
    # -------------------------------------------------------------------------
    # Agregamos un bloque para la variable 'inmunizado'
    index.append(("Nirsevimab Recipients", ""))
    tabla.append({})  # Fila vacía para separar visualmente

    # 1. Bloque categórico "Yes / No"
    for val_inm, label_inm in [(1, "Yes"), (0, "No")]:
        fila_inm = {}
        for g in grupos:
            subset = df[df['grupo_tabla'] == g]
            total = len(subset)
            if 'inmunizado' not in subset.columns:
                # Si la columna no existe, forzamos 0
                fila_inm[nombres_columnas[g]] = "- (No data)"
            else:
                cuenta = (subset['inmunizado'] == val_inm).sum()
                porcentaje = 100 * cuenta / total if total > 0 else 0
                fila_inm[nombres_columnas[g]] = f"{cuenta} ({porcentaje:.1f}%)"
        index.append(("Nirsevimab Recipients", label_inm))
        tabla.append(fila_inm)

    # 2. Bloque "Time immunized (days)"
    #    Resta 'fechaIng_any' - 'fecha_Inm'. Si no hay datos, "-", sino media en días

    # -------------------------------------------------------------------------
    # ---- Diagnósticos - 2 grupos: diag_labels y diag_agrupao
    # -------------------------------------------------------------------------
    diag_labels = {
        'A090': 'Infectious diarrhoea',
        'R681': 'Other general symptoms and signs',
        'R11X': 'Nausea and vomiting',
        'R104': 'Abdominal pain',
        'S099': 'Other head injuries',
        'R634': 'Feeding difficulties',
        'R633': 'Polydipsia',
        'N390': 'Urinary tract infection',
        'P599': 'Other neonatal jaundice',
    }

    diag_agrupao = {
        'A099': 'Other gastroenteritis and colitis',
        'A080': 'Other gastroenteritis and colitis',
        'A083': 'Other gastroenteritis and colitis',
        'A082': 'Other gastroenteritis and colitis',
        'A084': 'Other gastroenteritis and colitis',
        'A085': 'Other gastroenteritis and colitis',
        'A081': 'Other gastroenteritis and colitis',
    }

    diag_frecuencias = []
    # for diag_code, label in diag_labels.items():
    #     fila_diag = {}
    #     count_control = 0
    #     for g in grupos:
    #         if g == case_label:
    #             fila_diag[nombres_columnas[g]] = '-'  # Casos no exhiben "motivo"
    #         else:
    #             subset = df[(df['grupo_tabla'] == g) & (df['DIAG1'] == diag_code)]
    #             total = len(df[df['grupo_tabla'] == g])
    #             cuenta = len(subset)
    #             count_control = cuenta
    #             porcentaje = 100 * cuenta / total if total > 0 else 0
    #             fila_diag[nombres_columnas[g]] = f'{cuenta} ({porcentaje:.1f}%)'
    #     diag_frecuencias.append((count_control, ('Reason for hospital visit (controls)', label), fila_diag))

    # # Agrupación multiple: todos esos A09x a 'Other gastroenteritis and colitis'
    # # calculamos la suma total para los 'controls'
    # count_control_agrupao = 0
    # for diag_code, label in diag_agrupao.items():
    #     for g in grupos:
    #         if g == control_label:
    #             subset = df[(df['grupo_tabla'] == g) & (df['DIAG1'] == diag_code)]
    #             total = len(df[df['grupo_tabla'] == g])
    #             cuenta = len(subset)
    #             count_control_agrupao += cuenta

    # fila_agrupao = {}
    # porcentaje_agrupao = 100 * count_control_agrupao / total if total > 0 else 0
    # # Asignamos al control
    # fila_agrupao[nombres_columnas[control_label]] = f'{count_control_agrupao} ({porcentaje_agrupao:.1f}%)'
    # # Casos sin "motivo"
    # fila_agrupao[nombres_columnas[case_label]] = '-'
    # diag_frecuencias.append((
    #     count_control_agrupao,
    #     ('Reason for hospital visit (controls)', 'Other gastroenteritis and colitis'),
    #     fila_agrupao
    # ))

    # # Insertamos encabezado en la tabla
    # index.append(('Reason for hospital visit (controls)', ''))
    # tabla.append({})

    # # Ordenar por conteo y actualizar la tabla
    # for _, idx_, fila_ in sorted(diag_frecuencias, key=lambda x: x[0], reverse=True):
    #     index.append(idx_)
    #     tabla.append(fila_)

    multiindex = pd.MultiIndex.from_tuples(index, names=['Variable', 'Value'])
    return pd.DataFrame(tabla, index=multiindex)

def tabla_1_match_multi(
    df1=None, label1="base1",
    df2=None, label2="base2",
    df3=None, label3="base3",
    df4=None, label4="base4",
    df5=None, label5="base5",
    df6=None, label6="base6",
    id_col='RUN',
    group_col='Group',
    case_label='Caso',
    control_label='Control',
    decimales=2
):
    """
    Genera la tabla 1 para hasta 6 bases de datos. Si se pasa sólo df1, 
    reproduce la tabla de la función original. Para 2 a 6, concatena 
    horizontalmente las tablas.
    """
    # Lista de (df, label) dinámica
    df_list = []
    for df, lbl in [
        (df1, label1), (df2, label2), (df3, label3),
        (df4, label4), (df5, label5), (df6, label6)
    ]:
        if df is not None:
            df_list.append((df, lbl))

    if not df_list:
        print("No se recibió ninguna base de datos.")
        return None

    # Procesar cada base y renombrar columnas
    tables = []
    for df_x, lbl in df_list:
        single_table = _tabla_1_single(
            df_x,
            id_col=id_col,
            group_col=group_col,
            case_label=case_label,
            control_label=control_label,
            decimales=decimales
        )
        # MultiIndex para evitar colisiones
        single_table.columns = pd.MultiIndex.from_product([[lbl], single_table.columns])
        tables.append(single_table)

    # Devolver resultado según número de tablas
    if len(tables) == 1:
        return tables[0]
    else:
        return pd.concat(tables, axis=1)

    
# df1 = (pd.read_csv(path_data/'matched_data_VRS_general_francia.csv',index_col=0)
# .assign(time_inm = lambda x: (pd.to_datetime(x['fechaIng_any']) - pd.to_datetime(x['fechaInm'])).dt.days))
# df2 = (pd.read_csv(path_data/'matched_data_prematuros_all_born.csv',index_col=0)
# .assign(time_inm = lambda x: (pd.to_datetime(x['fechaIng_any']) - pd.to_datetime(x['fechaInm'])).dt.days))
# df3 = (pd.read_csv(path_data/'matched_data_cardiopats_1_all_born.csv',index_col=0)
# .assign(time_inm = lambda x: (pd.to_datetime(x['fechaIng_any']) - pd.to_datetime(x['fechaInm'])).dt.days))

df3 = (pd.read_csv(path_data/'matched_data_pali_all_born.csv',index_col=0)
       .assign(time_inm = lambda x: (pd.to_datetime(x['fechaIng_any']) - pd.to_datetime(x['fechaInm'])).dt.days))
df2 = (pd.read_csv(path_data/'matched_data_prematuros_no_pali_all_born.csv',index_col=0)
       .assign(time_inm = lambda x: (pd.to_datetime(x['fechaIng_any']) - pd.to_datetime(x['fechaInm'])).dt.days))
df1 = pd.concat([df2,df3])
df4 = (pd.read_csv(path_data/'matched_data_super_prematuros_all_born.csv',index_col=0)
       .assign(time_inm = lambda x: (pd.to_datetime(x['fechaIng_any']) - pd.to_datetime(x['fechaInm'])).dt.days))
df5 = (pd.read_csv(path_data/'matched_data_cardiopats_1_all_born.csv',index_col=0)
       .assign(time_inm = lambda x: (pd.to_datetime(x['fechaIng_any']) - pd.to_datetime(x['fechaInm'])).dt.days))


Tabla_1 = tabla_1_match_multi(
    df1=df1, label1="Pali + No-Pali",
    df2=df2, label2="Preterms-No-Pali",
    df3=df3, label3="Pali",
    df4=df4, label4="Super-Preterm",
    df5=df5, label5="Heart-congenital-disease",
    id_col='RUN',
    group_col='Group',
    case_label='Caso',
    control_label='Control',
    decimales=2
)

Tabla_1.to_csv(path_data/"tabla_1_match_multi_rev.csv")

with pd.ExcelWriter(path_data/"resultados_comparados_copy_nuevo.xlsx", engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    Tabla_1.to_excel(writer, sheet_name="Multi_tabla_1_v3_rev")
    
    
# import pandas as pd
# import numpy as np
# from scipy.stats import chi2_contingency

# def calcular_pvalores_tabla_1(tabla_1, grupo):
#     """
#     Calcula p-valores de comparación entre casos y controles dentro del grupo especificado.
#     """
#     caso_col = f"{grupo} (n = {tabla_1[grupo].columns[0].split('=')[1].strip(')')})"
#     control_col = f"{grupo} (n = {tabla_1[grupo].columns[1].split('=')[1].strip(')')})"

#     var_chi2 = {
#         "Age (days)": ["<90", "90-180", ">180"],
#         "Sex": ["Male", "Female"],
#         "Gestational age at birth (weeks)": ["24-28", "29-35"],
#         "Risk groups": ["Congenital heart disease", "Prematurity", "Prematurity Extreme"],
#         "Nirsevimab Recipients": ["Yes", "No"]
#     }

#     var_cont = "Birth weight (grams)"

#     # Extraer número desde string tipo "12 (34.5%)"
#     def extraer_n(text):
#         try:
#             return int(text.split('(')[0].strip())
#         except:
#             return np.nan

#     p_values = []

#     # Chi-cuadrado para categóricas
#     for var, categorias in var_chi2.items():
#         freqs = []
#         for cat in categorias:
#             fila = tabla_1.loc[(var, cat), (grupo, slice(None))]
#             n_caso = extraer_n(fila[(grupo, caso_col)])
#             n_control = extraer_n(fila[(grupo, control_col)])
#             freqs.append([n_caso, n_control])
#         try:
#             contingency = np.array(freqs)
#             _, p, _, _ = chi2_contingency(contingency.T)
#         except:
#             p = np.nan
#         p_values.append((var, p))

#     # Birth weight: no hay datos crudos, se deja como placeholder
#     p_values.append((var_cont, np.nan))

#     # Convertir a DataFrame
#     df_pvals = pd.DataFrame(p_values, columns=["Variable", "p-value"])
#     return df_pvals


# # Suponiendo que ya tienes `Tabla_1` generado
# df_pvals = calcular_pvalores_tabla_1(Tabla_1, grupo="Preterms-No-Pali")

# # Guardar en Excel
# with pd.ExcelWriter(path_data/"resultados_comparados_copy.xlsx", engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
#     df_pvals.to_excel(writer, sheet_name="tabla_p_values", index=False)