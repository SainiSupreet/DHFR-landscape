from _DFE import *
from _mutation_effect import *

path = '../DFE/'

# KS Test
background_bins = dfe_bins(list(fitness_data.values()),list(fitness_data.keys()),bin_size=0.05)

bin_values = list(background_bins.keys())
bin_values.sort()

bins = []
ks_p_mean = []
ks_p_median = []
min_ks = []
max_ks = []
samples = []

for binV in bin_values:
    print(f'Checking for bin {binV}')
    p_val = []

    count = 0
    for i1 in range(1,len(background_bins[binV])):
        for i2 in range(i1):

            count += 1
            if count >= 100000:
                break
            
            b1 = background_bins[binV][i1]
            e1 = []
            for mutation in generate_mutations(1):
                if mutation_feasibility(b1,mutation[0]) and get_mutant(b1,mutation[0])\
                in fitness_data.keys():
                    e1.append(fitness_data[get_mutant(b1,mutation[0])] - fitness_data[b1])
                    
            b2 = background_bins[binV][i2]
            e2 = []
            for mutation in generate_mutations(1):
                if mutation_feasibility(b2,mutation[0]) and get_mutant(b2,mutation[0])\
                in fitness_data.keys():
                    e2.append(fitness_data[get_mutant(b2,mutation[0])] - fitness_data[b2])
            
            p_val.append(stats.ks_2samp(e1,e2).pvalue)

        if count >= 100000:
            break
    if len(p_val) != 0:
        bins.append(binV)
        ks_p_mean.append(float(mean(p_val)))
        ks_p_median.append(float(median(p_val)))
        max_ks.append(float(max(p_val)))
        min_ks .append(float(min(p_val)))
        samples.append(count)

        df = pd.DataFrame({'p value':p_val})
        df.to_csv(f'../DFE/KS test/Background Fitness {binV}.csv',index = False)
    
ks_p_df = pd.DataFrame({'Bin':bins,
                        'mean p value':ks_p_mean,
                        'median p value':ks_p_median,
                        'Minimum p value':min_ks,
                        'Maximum p value':max_ks,
                        'Samples Examined':samples})
ks_p_df.to_csv(f'{path}KS test p values.csv',index = False)

