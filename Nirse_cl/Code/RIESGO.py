import sys
import os

parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

if parent_dir not in sys.path:
    sys.path.insert(0,parent_dir+'/src')

from config import config as cf
from config import pretreatment as pt

folder_path_input = cf.dataDirI
folder_path = cf.dataDirO
path_nacimientos = 'NAC_RNI_EGRESOS_ENTREGA_ISCI_11_04_2025_encr.csv'
output_text = '_04_2025'
nacimientos = pt.nacimientos(path_nacimientos,output_text)
path_nacimientos2 = f'nacimientos{output_text}.parquet'

path_egresos = 'egresos_5_total.parquet'
output_text2 = '_ensayo'
pt.riesgo(path_nacimientos2,path_egresos,output_text2)


