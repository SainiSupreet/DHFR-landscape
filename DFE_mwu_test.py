from _DFE import *
from _mutation_effect import *

path = '../DFE/'

# MWU Test
background_bins = dfe_bins(list(fitness_data.values()),list(fitness_data.keys()),bin_size=0.05)

bin_values = list(background_bins.keys())
bin_values.sort()

bins = []
mwu_mean = []
mwu_median = []
min_mwu = []
max_mwu = []
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
                    
            p_val.append(stats.mannwhitneyu(e1,e2).pvalue)

        if count >= 100000:
            break
    if len(p_val) != 0:
        bins.append(binV)
        mwu_mean.append(float(mean(p_val)))
        mwu_median.append(float(median(p_val)))
        max_mwu.append(float(max(p_val)))
        min_mwu.append(float(min(p_val)))
        samples.append(count)

        df = pd.DataFrame({'p value':p_val})
        df.to_csv(f'../DFE/MWU test/Background Fitness {binV}.csv',index = False)
    
mwu_df = pd.DataFrame({'Bin':bins,
                        'mean p value':mwu_mean,
                        'median p value':mwu_median,
                        'Minimum p value':min_mwu,
                        'Maximum p value':max_mwu,
                        'Samples Examined':samples})
mwu_df.to_csv(f'{path}MWU test p values.csv',index = False)

