from _DFE import *
from _mutation_effect import *
import scipy.io as sio

df = pd.read_csv('../DFE/Functional compiled data.csv') 
bins = dfe_bins(df['Background Fitness'],df['SC'],
                bin_size=0.05)
bin_values = list(bins.keys())
bin_values.sort()

data_array = []

x = bin_values
max_preset = max(df['SC'].values)
min_preset = min(df['SC'].values)

min_cut = 1
max_cut = 4

bin_sc = [round(float(min_preset+0.01*i),2) for i in range(0,int((max_preset - min_preset)//0.01 + 1))]
y = [float(i) for i in bin_sc]
z = []

for val in bin_sc:
    z.append([])
    for bins in bin_values:
        df = pd.read_csv(f'../DFE/Probability Distributions - functional/BF {bins}.csv')
        z[-1].append(float(df['PDF'][df['SC'] == val].values[0]))
    z[-1] = z[-1][min_cut:-max_cut]

x = x[min_cut:-max_cut]

matlab = '['

for l1 in z:
    for val in l1:
        matlab += str(val)
        matlab += ','
    matlab = matlab[:-1]
    matlab += ';'
matlab = matlab[:-1] + ']'

sio.savemat('../DFE/3d_data_functional.mat',dict(x=x,y=y,z=z))

df = pd.read_csv('../DFE/Non Functional compiled data.csv') 
bins = dfe_bins(df['Background Fitness'],df['SC'],
                bin_size=0.05)
bin_values = list(bins.keys())
bin_values.sort()

data_array = []

x = bin_values
max_preset = max(df['SC'].values)
min_preset = min(df['SC'].values)

min_cut = 1
max_cut = 1

bin_sc = [round(float(min_preset+0.01*i),2) for i in range(0,int((max_preset - min_preset)//0.01 + 1))]
y = [float(i) for i in bin_sc]
z = []

for val in bin_sc:
    z.append([])
    for bins in bin_values:
        df = pd.read_csv(f'../DFE/Probability Distributions - non functional/BF {bins}.csv')
        z[-1].append(float(df['PDF'][df['SC'] == val].values[0]))
    z[-1] = z[-1][min_cut:-max_cut]

x = x[min_cut:-max_cut]

matlab = '['

for l1 in z:
    for val in l1:
        matlab += str(val)
        matlab += ','
    matlab = matlab[:-1]
    matlab += ';'
matlab = matlab[:-1] + ']'

sio.savemat('../DFE/3d_data_non_functional.mat',dict(x=x,y=y,z=z))
