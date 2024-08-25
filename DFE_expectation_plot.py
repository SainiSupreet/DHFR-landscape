from bio_functions import *

plt.rcParams["figure.figsize"] = (50,50)
plt.rcParams['axes.linewidth'] = 20
plt.rcParams['lines.linewidth'] = 15
plt.rcParams['xtick.major.size']=50
plt.rcParams['xtick.major.width']=10
plt.rcParams['ytick.major.size']=50
plt.rcParams['ytick.major.width']=10

plt.rcParams['axes.labelpad'] = 60

plt.rcParams['boxplot.meanprops.linewidth'] = 10
plt.rcParams['boxplot.capprops.linewidth']= 10
plt.rcParams['boxplot.meanprops.markersize'] = 10
plt.rcParams['boxplot.whiskerprops.linewidth'] = 7
plt.rcParams['boxplot.flierprops.linewidth'] = 7
plt.rcParams['boxplot.boxprops.linewidth']=7

plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.constrained_layout.h_pad']=0.5
plt.rcParams['figure.constrained_layout.w_pad']=0.5
plt.rcParams['savefig.pad_inches'] = 0.5

plt.rcParams['grid.linewidth']=3

## All variants

n_bins = 10 # Number of bins in background fitness

max_bf = 1.395791684395265 # Maximum fitness value in data
min_bf = -1.17480288668069 # Minimum fitness value in data

bin_size = (max_bf - min_bf)/(n_bins)
bins = [(round(min_bf + i*bin_size,2) , round(min_bf + (i+1)*bin_size,2)) for i in range(n_bins)] #(min,max) of each bin

def find_bin(value,bin_size = bin_size,min_bf = min_bf):
    '''
    Returns the list index of the bin the value belongs to
    '''
    return int((value - min_bf)//bin_size)

background_bins = [[] for i in range(n_bins)]

## Segregate variants based on fitness

genes = list(fitness_data.keys())
fitness = list(fitness_data.values())

for i in range(len(genes)):
    background_bins[find_bin(fitness[i])].append(genes[i])

print(f'All variants data segregated:\n\
========================')
for i in range(n_bins):
    print(f'{bins[i]}:{len(background_bins[i])} variants')

## Checking DFEs
DFE_expectation_val = []
non_available_variants = set()

count = 0
for background_bin in background_bins:

    count += 1
    print(f'Checking DFE of bin: {bins[count-1]}')
    
    exp = []
    for variant in background_bin:

        bf = fitness_data[variant]
        sc = []

        for mutations in generate_mutations(1):

            mutation = mutations[0]

            if not mutation_feasibility(variant,mutation):
                continue

            mutant = get_mutant(variant,mutation)
            
            try:
                mf = fitness_data[mutant]
            except KeyError:
                non_available_variants.add(mutant)
                continue

            sc.append(mf-bf)

        exp.append(mean(sc))

    DFE_expectation_val.append(exp)
    print('Done')

plt.plot([0.5,n_bins + 0.5],[0,0],linewidth = 15,linestyle = '--',color = 'gray')

bp = plt.boxplot(DFE_expectation_val,showmeans = True,meanline = True,
                 showfliers = False,
                 medianprops=dict(linestyle='-', linewidth=15,
                                  color='black'),
                 meanprops=dict(linestyle='--', linewidth=15,
                                color='red'))

plt.legend([bp['means'][0],bp['medians'][0]],
           ['Mean','Median'],frameon=False,fontsize=150)

plt.xticks(range(1,n_bins+1),[round(mean([bbin[0],bbin[1]]),2) for bbin in bins],
           fontsize = 150)

plt.yticks(fontsize = 150)
plt.locator_params(axis='both',nbins=6)

plt.xlabel('Background fitness',fontsize = 175)
plt.ylabel('Selection coefficient',fontsize = 175)
plt.savefig(f'../Plots/DFE/Mean DFE (All variants)')

plt.close()

## Functional variants

n_bins = 10 # Number of bins in background fitness

max_bf = 1.395791684395265 # Maximum fitness value in data
min_bf = threshold # Minimum fitness value in data

bin_size = (max_bf - min_bf)/(n_bins)
bins = [(round(min_bf + i*bin_size,2) , round(min_bf + (i+1)*bin_size,2)) for i in range(n_bins)] #(min,max) of each bin

def find_bin(value,bin_size = bin_size,min_bf = min_bf):
    '''
    Returns the list index of the bin the value belongs to
    '''
    return int((value - min_bf)//bin_size)

background_bins = [[] for i in range(n_bins)]

## Segregate variants based on fitness

genes = list(fitness_data.keys())
fitness = list(fitness_data.values())

for i in range(len(genes)):

    # Remove non functional variants
    if fitness[i] < threshold:
        continue
    
    background_bins[find_bin(fitness[i])].append(genes[i])

print(f'Functional variants data segregated:\n\
========================')
for i in range(n_bins):
    print(f'{bins[i]}:{len(background_bins[i])} variants')

## Checking DFEs
DFE_expectation_val = []
non_available_variants = set()

count = 0
for background_bin in background_bins:

    count += 1
    print(f'Checking DFE of bin: {bins[count-1]}')
    
    exp = []
    for variant in background_bin:

        bf = fitness_data[variant]
        sc = []

        for mutations in generate_mutations(1):

            mutation = mutations[0]

            if not mutation_feasibility(variant,mutation):
                continue

            mutant = get_mutant(variant,mutation)
            
            try:
                mf = fitness_data[mutant]
            except KeyError:
                non_available_variants.add(mutant)
                continue

            sc.append(mf-bf)

        exp.append(mean(sc))

    DFE_expectation_val.append(exp)
    print('Done')

plt.plot([0.5,n_bins + 0.5],[0,0],linewidth = 15,linestyle = '--',color = 'gray')

bp = plt.boxplot(DFE_expectation_val,showmeans = True,meanline = True,
                 showfliers = False,
                 medianprops=dict(linestyle='-', linewidth=15,
                                  color='black'),
                 meanprops=dict(linestyle='--', linewidth=15,
                                color='red'))

plt.legend([bp['means'][0],bp['medians'][0]],
           ['Mean','Median'],frameon=False,fontsize=150)

plt.xticks(range(1,n_bins+1),[round(mean([bbin[0],bbin[1]]),2) for bbin in bins],
           fontsize = 150)

plt.yticks(fontsize = 150)
plt.locator_params(axis='both',nbins=6)

plt.xlabel('Background fitness',fontsize = 175)
plt.ylabel('Selection coefficient',fontsize = 175)
plt.savefig(f'../Plots/DFE/Mean DFE (Functional variants)')

plt.close()

## Non functional variants

n_bins = 10 # Number of bins in background fitness

max_bf = threshold # Maximum fitness value in data
min_bf = -1.17480288668069 # Minimum fitness value in data

bin_size = (max_bf - min_bf)/(n_bins)
bins = [(round(min_bf + i*bin_size,2) , round(min_bf + (i+1)*bin_size,2)) for i in range(n_bins)] #(min,max) of each bin

def find_bin(value,bin_size = bin_size,min_bf = min_bf):
    '''
    Returns the list index of the bin the value belongs to
    '''
    return int((value - min_bf)//bin_size)

background_bins = [[] for i in range(n_bins)]

## Segregate variants based on fitness

genes = list(fitness_data.keys())
fitness = list(fitness_data.values())

for i in range(len(genes)):

    if fitness[i] > threshold:
        continue
    
    background_bins[find_bin(fitness[i])].append(genes[i])

print(f'Non functional variants data segregated:\n\
========================')
for i in range(n_bins):
    print(f'{bins[i]}:{len(background_bins[i])} variants')

## Checking DFEs
DFE_expectation_val = []
non_available_variants = set()

count = 0
for background_bin in background_bins:

    count += 1
    print(f'Checking DFE of bin: {bins[count-1]}')
    
    exp = []
    for variant in background_bin:

        bf = fitness_data[variant]
        sc = []

        for mutations in generate_mutations(1):

            mutation = mutations[0]

            if not mutation_feasibility(variant,mutation):
                continue

            mutant = get_mutant(variant,mutation)
            
            try:
                mf = fitness_data[mutant]
            except KeyError:
                non_available_variants.add(mutant)
                continue

            sc.append(mf-bf)

        exp.append(mean(sc))

    DFE_expectation_val.append(exp)
    print('Done')

plt.plot([0.5,n_bins + 0.5],[0,0],linewidth = 15,linestyle = '--',color = 'gray')

bp = plt.boxplot(DFE_expectation_val,showmeans = True,meanline = True,
                 showfliers = False,
                 medianprops=dict(linestyle='-', linewidth=15,
                                  color='black'),
                 meanprops=dict(linestyle='--', linewidth=15,
                                color='red'))

plt.legend([bp['means'][0],bp['medians'][0]],
           ['Mean','Median'],frameon=False,fontsize=150)

plt.xticks(range(1,n_bins+1),[round(mean([bbin[0],bbin[1]]),2) for bbin in bins],
           fontsize = 150)

plt.yticks(fontsize = 150)
plt.locator_params(axis='both',nbins=6)

plt.xlabel('Background fitness',fontsize = 175)
plt.ylabel('Selection coefficient',fontsize = 175)
plt.savefig(f'../Plots/DFE/Mean DFE (Non functional variants)')

plt.close()
