import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import sys
import os
from dateutil.relativedelta import relativedelta

parent_dir = os.getcwd()
if parent_dir not in sys.path: sys.path.insert(0,parent_dir+'/src')

from config import config as cf
from config import functions as f


folder_path_input = cf.dataDirI2
folder_path = cf.dataDirO2
egresos = pd.read_parquet(folder_path+'/egresos_all.parquet')
egresos_2024 = pd.read_parquet(folder_path+'/egresos_all_2024.parquet')
df5 = pd.read_parquet(folder_path+'/egresos_5.parquet')
df5_2024 = pd.read_parquet(folder_path+'/egresos_5_2024.parquet')
dfrc = pd.read_parquet(f'{folder_path}/egresos_respiratory_causes.parquet')
dfrc_2024 = pd.read_parquet(f'{folder_path}/egresos_respiratory_causes_2024.parquet')
pd.concat([egresos,egresos_2024]).to_parquet(f'{folder_path}/egresos_total.parquet', engine='pyarrow', compression='snappy')
pd.concat([df5,df5_2024]).to_parquet(f'{folder_path}/egresos_5_total.parquet', engine='pyarrow', compression='snappy')
pd.concat([dfrc,dfrc_2024]).to_parquet(f'{folder_path}/egresos_respiratory_causes_total.parquet', engine='pyarrow', compression='snappy')



