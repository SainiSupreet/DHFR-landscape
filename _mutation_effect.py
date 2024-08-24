from bio_functions import *

def study_mutation(mutation_list:list,background_status = ['non functional','functional','all'],
                   consider_threshold = False,save_data = False,path = './',
                   fitness_data = fitness_data,threshold = threshold):
    '''
    Study the effect of mutation in all backgrounds
    '''
    if background_status not in ['non functional','functional','all']:
        raise ValueError('Invalid Background demand')
    
    if len(mutation_list)!=1:
        raise ValueError('Exactly 1 mutation is required')
    
    mutation = mutation_list[0]
    
    compiled_data = {
        'Background':[],
        'Background Fitness':[],
        'Mutant':[],
        'Mutant Fitness':[],
        'SC':[],
    }
    
    for background in fitness_data.keys():

        if (fitness_data[background] < threshold and background_status == 'functional') or \
           (fitness_data[background] > threshold and background_status == 'non functional'):
            continue
        
        if not mutation_feasibility(background,mutation):
            continue
        
        bf = fitness_data[background]
        mutant = get_mutant(background,mutation)
        
        try:
            mf = fitness_data[mutant]
        except KeyError:
            continue
        
        if consider_threshold and (fitness_data[background]<=threshold and fitness_data[mutant]<=threshold):
            continue
        
        s = mf - bf
        
        compiled_data['Background'].append(background)
        compiled_data['Background Fitness'].append(bf)
        compiled_data['Mutant'].append(mutant)
        compiled_data['Mutant Fitness'].append(mf)
        compiled_data['SC'].append(s)
    
    if save_data:
        
        compiled_df = pd.DataFrame(data = compiled_data)
        compiled_df.to_csv(f'{path}{mutation[0]} at {mutation[1]+1} mutation data.csv',
                           index = False)
    
    return compiled_data

# Expanded list of all mutations - including separate cases for base
mutations = generate_mutations(1)
expanded_mutations = []

for mut in mutations:
    mutation = mut[0]
    for base in nucleotides:
        if mutation[0]!=base:
            expanded_mutations.append((base,mutation[0],mutation[1]))
