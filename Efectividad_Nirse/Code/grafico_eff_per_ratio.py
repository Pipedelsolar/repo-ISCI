import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import re
import numpy as np
from pathlib import Path

# 1) Leer el archivo Excel (todas las hojas)
path_actual = Path.cwd()
path_data = path_actual.parent / 'Data'
file_path = path_data / "results_compare_python.xlsx"
excel = pd.ExcelFile(file_path)
sheet_names = excel.sheet_names

# 2) Función para parsear la columna "efectividad"
def parse_efectividad(s):
    """
    Retorna (mean, lw, up) como floats
    Ejemplo: "75.65 (58.1; 85.85)"
    """
    pattern = r"([\d\.]+)\s*\(([\d\.]+);\s*([\d\.]+)\)"
    match = re.search(pattern, s)
    if match:
        mean_val = float(match.group(1))
        lw_val = float(match.group(2))
        up_val = float(match.group(3))
        return mean_val, lw_val, up_val
    else:
        return None, None, None

# 3) Recolectar datos en un DataFrame maestro
master_rows = []

for sheet in sheet_names:
    df_sheet = excel.parse(sheet_name=sheet)
    # Asegurarte de que la columna se llame "efectividad" y también "filtro"
    if "efectividad" not in df_sheet.columns or "filtro" not in df_sheet.columns:
        continue

    for i in range(len(df_sheet)):
        val = df_sheet.loc[i, "efectividad"]
        mean, lw, up = parse_efectividad(str(val))
        master_rows.append({
            "sheet": sheet,
            "row": i,
            "mean": mean,
            "lw": lw,
            "up": up,
            # Guarda también el filtro
            "filtro": df_sheet.loc[i, "filtro"]
        })

df_master = pd.DataFrame(master_rows)

# 4) Graficar y guardar
unique_rows = df_master['row'].unique()
for row_idx in unique_rows:
    sub = df_master[df_master["row"] == row_idx].copy()
    #df_vrs_match_casesub = sub.sort_values("sheet")

    x = range(len(sub))
    means = sub["mean"]
    y_lower = means - sub["lw"]
    y_upper = sub["up"] - means

    if not sub.empty:
        filtro_name = sub.iloc[0]["filtro"]
    else:
        filtro_name = f"Fila_{row_idx}"

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar(
        x,
        means,
        yerr=[y_lower, y_upper],
        fmt='o',
        ecolor='blue',
        capsize=5,
        label="IC 95%"
    )
    ax.set_title(f"{filtro_name} - Efectividad con IC")
    ax.set_xlabel("Nombre de la hoja")
    ax.set_ylabel("Efectividad (%)")

    ax.set_xticks(x)
    ax.set_xticklabels(sub["sheet"], rotation=45, ha='right')

    plt.tight_layout()

    # Guardar la figura en un archivo PNG (o PDF, JPG, etc.)
    # Ejemplo: "filtro_name_EFIC.png" (reemplaza espacios/caracteres especiales si es necesario)
    safe_filtro_name = str(filtro_name).replace(" ", "_")  # quitar espacios
    save_path = path_actual / f"{safe_filtro_name}_efect.png"
    plt.savefig(save_path, dpi=300)  # alta resolución

    # Si quieres ver la figura en pantalla:
    plt.show()
    # O si no quieres mostrarla, podrías hacer plt.close(fig)



















# import pandas as pd
# import matplotlib.pyplot as plt
# import re
# import numpy as np
# from pathlib import Path

# # ------------------------------------------------------------------------------
# # 1) Leer el archivo Excel (todas las hojas)
# # ------------------------------------------------------------------------------
# path_actual = Path.cwd()
# path_data = path_actual.parent / 'Data'
# file_path = path_data / "results_compare_python.xlsx"
# excel = pd.ExcelFile(file_path)
# sheet_names = excel.sheet_names

# # ------------------------------------------------------------------------------
# # 2) Función para parsear la columna "efectividad", tipo "75.65 (58.1; 85.85)"
# # ------------------------------------------------------------------------------
# def parse_efectividad(s):
#     pattern = r"([\d\.]+)\s*\(([\d\.]+);\s*([\d\.]+)\)"
#     match = re.search(pattern, s)
#     if match:
#         mean_val = float(match.group(1))
#         lw_val = float(match.group(2))
#         up_val = float(match.group(3))
#         return mean_val, lw_val, up_val
#     else:
#         return None, None, None  # Problema si no coincide el formato

# # ------------------------------------------------------------------------------
# # 3) Construir DataFrame maestro (df_master) con filas de cada hoja
# # ------------------------------------------------------------------------------
# master_rows = []

# for sheet in sheet_names:
#     df_sheet = excel.parse(sheet_name=sheet)

#     # Asegúrate que "efectividad" y "filtro" estén en la hoja
#     if "efectividad" not in df_sheet.columns or "filtro" not in df_sheet.columns:
#         continue

#     for i in range(len(df_sheet)):
#         val = df_sheet.loc[i, "efectividad"]
#         mean, lw, up = parse_efectividad(str(val))
#         master_rows.append({
#             "sheet": sheet,
#             "row": i,   # Subgrupo (fila)
#             "mean": mean,
#             "lw": lw,
#             "up": up,
#             "filtro": df_sheet.loc[i, "filtro"]
#         })

# df_master = pd.DataFrame(master_rows)

# # ------------------------------------------------------------------------------
# # 4) Graficar todas las filas en un solo gráfico
# # ------------------------------------------------------------------------------
# # Listado único de hojas (ordenadas) y subgrupos (filas)
# sheet_list = sorted(df_master["sheet"].unique())
# unique_rows = sorted(df_master["row"].unique())

# # Creamos la figura
# fig, ax = plt.subplots(figsize=(12, 6))

# # Factor para aumentar el espacio entre hojas
# spacing_factor = 3.0  

# for idx, row_idx in enumerate(unique_rows):
#     # Filtramos solo la fila row_idx
#     sub = df_master[df_master["row"] == row_idx].copy()

#     # Convertir la hoja a una posición en el eje X
#     sub["x"] = sub["sheet"].apply(lambda s: sheet_list.index(s) * spacing_factor)
#     sub = sub.sort_values("x")

#     x_vals = sub["x"].values
#     means = sub["mean"].values
#     y_lower = means - sub["lw"].values  # Distancia entre media e IC inferior
#     y_upper = sub["up"].values - means  # Distancia entre media e IC superior

#     # Desplazamiento adicional para no superponer con otras filas
#     offset = idx * 0.4
#     x_offset = x_vals + offset

#     # Leyenda: usamos el valor de la primera fila (mismo 'filtro' en todas las hojas)
#     if not sub.empty:
#         legend_label = sub.iloc[0]["filtro"]
#     else:
#         legend_label = f"Fila {row_idx}"

#     # Graficar con errorbar
#     ax.errorbar(
#         x_offset,
#         means,
#         yerr=[y_lower, y_upper],
#         fmt='o',
#         capsize=5,
#         label=legend_label
#     )

# # Título y etiquetas de ejes
# ax.set_title("Efectividad con IC - Todas las filas en un gráfico (Mayor Separación)")
# ax.set_xlabel("Hoja de Excel")
# ax.set_ylabel("Efectividad (%)")

# # Ajustar las ubicaciones de ticks en el eje X
# x_positions = [i * spacing_factor for i in range(len(sheet_list))]
# ax.set_xticks(x_positions)
# ax.set_xticklabels(sheet_list, rotation=45, ha='right')

# ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
# plt.tight_layout()
# plt.show()

