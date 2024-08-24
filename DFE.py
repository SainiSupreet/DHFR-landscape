from _DFE import *
from _mutation_effect import *

all_data = {
    'Background Fitness':[],
    'Mutant Fitness':[],
    'SC':[]
    }
fun_data = {
    'Background Fitness':[],
    'Mutant Fitness':[],
    'SC':[]
    }
nf_data = {
    'Background Fitness':[],
    'Mutant Fitness':[],
    'SC':[]
    }

total_mutations = len(generate_mutations(1))
count = 0

path = '../DFE/'

for mutation in generate_mutations(1):
    
    data = study_mutation(mutation,
                          background_status = 'all',
                          save_data = False,
                          consider_threshold = False)
    fun = study_mutation(mutation,
                          background_status = 'functional',
                          save_data = False,
                          consider_threshold = False)
    nf = study_mutation(mutation,
                          background_status = 'non functional',
                          save_data = False,
                          consider_threshold = False)

    for key in all_data.keys():
        
        all_data[key]+=data[key]
        fun_data[key]+=fun[key]
        nf_data[key]+=nf[key]
        
    count += 1
    if count in [total_mutations//(100/i) for i in range(10,101,10)]:
        print(f'{round(100*count/total_mutations,0)}% mutations checked')

all_df = pd.DataFrame(all_data)
all_df = all_df.sort_values(by = 'Background Fitness')
all_df.to_csv(f'{path}All compiled data.csv',index = False)

fun_df = pd.DataFrame(fun_data)
fun_df = fun_df.sort_values(by = 'Background Fitness')
fun_df.to_csv(f'{path}Functional compiled data.csv',index = False)

nf_df = pd.DataFrame(nf_data)
nf_df = nf_df.sort_values(by = 'Background Fitness')
nf_df.to_csv(f'{path}Non functional compiled data.csv',index = False)

bin_values = dfe_bins(fun_df['SC'],fun_df['SC'],bin_size = 0.01)
bin_values = list(bin_values.keys())
bin_values.sort()

all_df = pd.read_csv('../DFE/All compiled data.csv')
fun_df = pd.read_csv('../DFE/Functional compiled data.csv')
nf_df = pd.read_csv('../DFE/Non Functional compiled data.csv')

# Functional Probability Distributions
bins = dfe_bins(fun_df['Background Fitness'],fun_df['SC'],
                bin_size=0.05)
bin_values = list(bins.keys())
bin_values.sort()

for val in bin_values:
    data = bin_data(bins[val],bin_size=0.01,
                    max_preset = max(fun_df['SC'].values),
                    min_preset = min(fun_df['SC'].values))
    norm_coeff = sum(list(data.values()))*0.01
    prob_dist = {i:data[i]/norm_coeff for i in data.keys()}
    
    dfe_df = pd.DataFrame(data = {'SC':prob_dist.keys(),'PDF':prob_dist.values()}).sort_values(by = 'SC')
    dfe_df.to_csv(f'{path}Probability Distributions - functional/BF {val}.csv',
                  index = False)

# Non Functional Probability Distributions
bins = dfe_bins(nf_df['Background Fitness'],nf_df['SC'],
                bin_size=0.05)
bin_values = list(bins.keys())
bin_values.sort()

for val in bin_values:
    data = bin_data(bins[val],bin_size=0.01,
                    max_preset = max(nf_df['SC'].values),
                    min_preset = min(nf_df['SC'].values))
    norm_coeff = sum(list(data.values()))*0.01
    prob_dist = {i:data[i]/norm_coeff for i in data.keys()}
    
    dfe_df = pd.DataFrame(data = {'SC':prob_dist.keys(),'PDF':prob_dist.values()}).sort_values(by = 'SC')
    dfe_df.to_csv(f'{path}Probability Distributions - non functional/BF {val}.csv',
                  index = False)

# All Probability Distributions
bins = dfe_bins(all_df['Background Fitness'],all_df['SC'],
                bin_size=0.05)
bin_values = list(bins.keys())
bin_values.sort()

for val in bin_values:
    data = bin_data(bins[val],bin_size=0.01,
                    max_preset = max(all_df['SC'].values),
                    min_preset = min(all_df['SC'].values))
    norm_coeff = sum(list(data.values()))*0.01
    prob_dist = {i:data[i]/norm_coeff for i in data.keys()}
    
    dfe_df = pd.DataFrame(data = {'SC':prob_dist.keys(),'PDF':prob_dist.values()}).sort_values(by = 'SC')
    dfe_df.to_csv(f'{path}Probability Distributions - all/BF {val}.csv',
                  index = False)
