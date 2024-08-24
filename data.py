from libraries import *

file = 'fitness_data_wt.rds'
result = pyreadr.read_r(file)
df1 = result[None]
variants = df1.SV
fitnesses = df1.m

fitness_data = {}

for i in range(len(variants)):
    fitness_data[variants[i]] = fitnesses[i]

print(len(fitness_data),'data compiled')

nucleotides = ['A','T','G','C']

threshold = -0.507774

no_epi_margin = 0.05
