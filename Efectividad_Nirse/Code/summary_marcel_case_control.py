import pandas as pd
import numpy as np

#results['max_cases']=df_case.shape[0]

# Lista para almacenar los DataFrames de cada experimento con el nombre del experimento
experiment_dfs = []
experiment_dfs_marcel = []

# Crear una función para evaluar el criterio
def check_signs(col1, col2):
    signs = np.sign([col1, col2])  # Obtiene los signos de ambos valores
    if all(signs > 0) or all(signs < 0):
        return 1  # Ambos positivos o ambos negativos
    else:
        return 0  # Signos intercalados



for result in list_experiments:
    # Obtener el DataFrame de resultados
    results_df_transposed =(
        result['results_df'].T
        .assign(
            str_values = lambda df : f'{str(df["Efectividad"].values[0].round(2))} ({str(df["0.025 Conf Interval"].values[0].round(2))}; {str(df["0.975 Conf Interval"].values[0].round(2))})',
            n_matched_cases = result['df_matched'].Group.nunique(),
            porc_n_matched_cases = result['df_matched'].Group.nunique()/result['max_cases'],
            pseudo_p_value = lambda df : check_signs(df['0.025 Conf Interval'], df['0.975 Conf Interval'])
        )
        .assign(
                method = result['method'],
                subset = result['subset'],
                region = result['region']
            )
        .pipe(lambda df:df.set_index(['region','method','subset']))
        )
    # Añadir el DataFrame a la lista
    experiment_dfs.append(results_df_transposed)

# Concatenar todos los DataFrames en uno solo
all_results_df = pd.concat(experiment_dfs, ignore_index=False)

# Mostrar el DataFrame concatenado
display(all_results_df)




marcel_summary= (all_results_df[['Efectividad']]
#.reset_index()
.unstack(-2)
.unstack(-1)

)

display(marcel_summary)

def generate_marcel_summary(all_results_df, porc_threshold=0.1, n_cases_threshold=5):
    """
    Genera marcel_summary filtrando por:
      - porc_n_matched_cases menor a un umbral (se reemplazan los valores con NaN).
      - n_matched_cases menor al umbral (se reemplazan los valores con NaN).
      - pseudo_p_value distinto de 1 (se reemplazan los valores con NaN).
    
    Parámetros:
    - all_results_df (DataFrame): DataFrame principal con resultados de experimentos.
    - porc_threshold (float): Umbral para filtrar porc_n_matched_cases.
    - n_cases_threshold (int): Umbral para filtrar n_matched_cases.

    Retorna:
    - DataFrame: Resumen con las condiciones aplicadas.
    """
    # Crear copia del DataFrame para no modificar el original
    filtered_df = all_results_df.copy()

    # Asignar NaN donde porc_n_matched_cases < porc_threshold
    filtered_df.loc[filtered_df['porc_n_matched_cases'] < porc_threshold, 'Efectividad'] = np.nan

    # Asignar NaN donde n_matched_cases < n_cases_threshold
    filtered_df.loc[filtered_df['n_matched_cases'] < n_cases_threshold, 'Efectividad'] = np.nan

    # Asignar NaN donde pseudo_p_value != 1
    filtered_df.loc[filtered_df['pseudo_p_value'] != 1, 'Efectividad'] = np.nan

    # Generar el resumen marcel_summary
    marcel_summary = (
        filtered_df[['Efectividad']]
        .unstack(-2)  # Desapilar por 'method'
        .unstack(-1)  # Desapilar por 'subset'
    )

    return marcel_summary


# Generar el resumen con las condiciones dadas
porc_threshold = 0.1  # Umbral para porc_n_matched_cases
n_cases_threshold = 5  # Umbral para n_matched_cases

marcel_summary_filtered = generate_marcel_summary(
    all_results_df, 
    porc_threshold=porc_threshold, 
    n_cases_threshold=n_cases_threshold
)

# Mostrar el resultado

porc_n_matched_cases= (all_results_df[['porc_n_matched_cases']]
#.reset_index()
.unstack(-2)
.unstack(-1)

)

n_matched_cases= (all_results_df[['n_matched_cases']]
#.reset_index()
.unstack(-2)
.unstack(-1)

)
display(marcel_summary_filtered)
display(porc_n_matched_cases)
display(n_matched_cases)

#graficar
import seaborn as sns
import matplotlib.pyplot as plt

# Graficar el heatmap del DataFrame numérico
plt.figure(figsize=(14, 8))
sns.heatmap(
    marcel_summary_filtered, 
    annot=True, 
    fmt=".2f", 
    cmap="RdBu_r",  # Azul para valores altos, rojo para valores bajos
    cbar=True, 
    linewidths=0.5,
    center=0  # Centrar la escala en 0
)
plt.title("Heatmap of Numeric DataFrame (High: Blue, Low: Red)")
plt.xlabel("Methods and Subsets")
plt.ylabel("Regions")
plt.tight_layout()
plt.show()

# Graficar el heatmap del DataFrame binario
plt.figure(figsize=(14, 8))
sns.heatmap(
    binary_df, 
    annot=True, 
    fmt="d", 
    cmap="Blues",  # Azul para activaciones (1), blanco para inactivaciones (0)
    cbar=False, 
    linewidths=0.5
)
plt.title("Binary Activation Heatmap")
plt.xlabel("Methods and Subsets")
plt.ylabel("Regions")
plt.tight_layout()
plt.show()