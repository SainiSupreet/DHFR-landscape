from _mutation_effect import *

for mutation in generate_mutations(1):
    path = '../Mutation effect - all points/'
    study_mutation(mutation,background_status = 'non functional',consider_threshold = False,save_data = True,path = path+'non functional/')
    study_mutation(mutation,background_status = 'functional',consider_threshold = False,save_data = True,path = path+'functional/')
    study_mutation(mutation,background_status = 'all',consider_threshold = False,save_data = True,path = path+'all/')
    
    path = '../Mutation effect - functional/'
    study_mutation(mutation,background_status = 'non functional',consider_threshold = True,save_data = True,path = path+'non functional/')
    study_mutation(mutation,background_status = 'functional',consider_threshold = True,save_data = True,path = path+'functional/')
    study_mutation(mutation,background_status = 'all',consider_threshold = True,save_data = True,path = path+'all/')

# Studying expanded mutations - for linear regression
count = 0
for mutation in expanded_mutations:

    data = {'Background':[],
            'Mutant':[],
            'Background Fitness':[],
            'SC':[]}
    
    for variant in fitness_data.keys():
        if variant[mutation[2]] != mutation[0]:
            continue

        mutant = get_mutant(variant,mutation[1:])

        try:
            mf = fitness_data[mutant]
            bf = fitness_data[variant]
            sc = mf - bf
        except KeyError:
            continue

        data['Background'].append(variant)
        data['Mutant'].append(mutant)
        data['Background Fitness'].append(bf)
        data['SC'].append(sc)

    df = pd.DataFrame(data)
    df.to_csv(f'../Mutation effect - all points/Expanded mutations - all (for linear regression)/{mutation[0]}->{mutation[1]} at {mutation[2]+1}.csv',index = False)

    count += 1
    if count in [108//(100/i) for i in range(10,101,10)]:
        print(f'{round(100*count/108,0)}% mutations checked')

folders = ['../Mutation effect - all points/']

for folder in folders:
    path = folder+'Expanded mutations - all (for linear regression)/'

    data_fit = {'Mutation base':[],
                'Mutation locus':[],
                'Linear slope':[],
                'Intercept':[],
                'Pivot point (background fitness)':[],
                'score (R^2)':[]
                }

    for mutation in expanded_mutations:

        background = mutation[0]
        mut = mutation[1]
        locus = mutation[2]+1

        file_name = f'{background}->{mut} at {locus}'

        data_fit['Mutation base'].append(f'{background} -> {mut}')
        data_fit['Mutation locus'].append(locus)

        #Fitting a linear model for bf vs. sc
        data = pd.read_csv(f'{path}{file_name}.csv')
        bf = data['Background Fitness'].values.reshape([len(data['Background Fitness']),1])
        sc = data['SC'].values.reshape([len(data['Background Fitness']),1])
        x = LR().fit(bf,sc)

        slope = x.coef_[0][0]
        intercept = x.intercept_[0]
        score = x.score(bf,sc)

        data_fit['Linear slope'].append(slope)
        data_fit['Intercept'].append(intercept)
        data_fit['Pivot point (background fitness)'].append(-intercept/slope)
        data_fit['score (R^2)'].append(score)

    # Checkpoint
    df = pd.DataFrame(data_fit)

    nf = {'BF':[],'SC':[]}
    fun = {'BF':[],'SC':[]}

    # Compiled effect of all mutations
    for mutations in generate_mutations(1):
        mutation = mutations[0]

        fun_data = pd.read_csv(f'{folder}functional/{mutation[0]} at {mutation[1]+1} mutation data.csv')
        nf_data = pd.read_csv(f'{folder}non functional/{mutation[0]} at {mutation[1]+1} mutation data.csv')

        fun['BF']+=list(fun_data['Background Fitness'].values)
        fun['SC']+=list(fun_data['SC'].values)

        nf['BF']+=list(nf_data['Background Fitness'].values)
        nf['SC']+=list(nf_data['SC'].values)

    fun_bf = np.array(fun['BF']).reshape([len(fun['SC']),1])
    fun_sc = np.array(fun['SC']).reshape([len(fun['SC']),1])

    nf_bf = np.array(nf['BF']).reshape([len(nf['SC']),1])
    nf_sc = np.array(nf['SC']).reshape([len(nf['SC']),1])

    fun_fit = LR().fit(fun_bf,fun_sc)
    nf_fit = LR().fit(nf_bf,nf_sc)

    # Functional variants
    slope = fun_fit.coef_[0][0]
    intercept = fun_fit.intercept_[0]
    score = fun_fit.score(fun_bf,fun_sc)
    
    data_fit['Mutation base'].append('Functional')
    data_fit['Mutation locus'].append('All')
    data_fit['Linear slope'].append(slope)
    data_fit['Intercept'].append(intercept)
    data_fit['Pivot point (background fitness)'].append(-intercept/slope)
    data_fit['score (R^2)'].append(score)

    # Non functional variants
    slope = nf_fit.coef_[0][0]
    intercept = nf_fit.intercept_[0]
    score = nf_fit.score(nf_bf,nf_sc)
    
    data_fit['Mutation base'].append('Non functional')
    data_fit['Mutation locus'].append('All')
    data_fit['Linear slope'].append(slope)
    data_fit['Intercept'].append(intercept)
    data_fit['Pivot point (background fitness)'].append(-intercept/slope)
    data_fit['score (R^2)'].append(score)

    df2 = pd.DataFrame(data_fit)
    df2.to_csv(f'{folder}Linear regression.csv',index = False)

    compiled_data = {'Mean pivot point':np.mean(df['Pivot point (background fitness)']),
                     'stdev of pivot point':np.std(df['Pivot point (background fitness)']),
                     'Cutoff fitness for functional variants':threshold,
                     'Mean R^2':np.mean(df['score (R^2)']),
                     'stdev of R^2':np.std(df['score (R^2)']),
                     }

    for i in expanded_mutations:
        compiled_data[f'R^2 of {i[0]} -> {i[1]} at {i[2]+1}']=df[df['Mutation base']==f'{i[0]} -> {i[1]}']\
                                                    [df['Mutation locus']==i[2]+1]\
                                                    ['score (R^2)'].values[0]

    compiled_data[f'R^2 of all high fitness backgrounds'] = fun_fit.score(fun_bf,fun_sc)
    compiled_data[f'R^2 of all low fitness backgrounds'] = nf_fit.score(nf_bf,nf_sc)
    
    compiled_df = pd.DataFrame({'Parameter':compiled_data.keys(),
                                'Value':compiled_data.values()})
    compiled_df.to_csv(f'{folder}Compiled results of linear regression.csv',index = False)
