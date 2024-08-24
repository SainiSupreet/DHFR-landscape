from bio_functions import *

data = pd.read_csv('../Beta lactamase/b_lactamase_mic.csv')

fitness = {i:data[data['Genotype'] == i]['MIC'].values[0] for i in data['Genotype'].values}
mutations = [('C',i) for i in range(5)]

def epistasis(background:str,mut1,mut2,fitness = fitness):
    '''
    Returns nature of epistasis in the given background
    '''
    if mut1 == mut2:
        raise ValueError('both mutations are the same')
    
    m1 = get_mutant(background,mut1)
    m2 = get_mutant(background,mut2)
    m12 = get_mutant(get_mutant(background,mut1),mut2)

    bf = fitness[background]
    f1 = fitness[m1]
    f2 = fitness[m2]
    f12 = fitness[m12]

    s1 = f1 - bf
    s2 = f2 - bf
    s12 = f12 - bf

    if s12*(s1 + s2) < 0:
        return 'Sign'
    elif s12 > (s1 + s2):
        return 'Positive'
    elif s12 < (s1 + s2):
        return 'Negative'
    else:
        return 'No'

def epistasis_change(background,mut1,mut2,mutations = mutations):
    ans = {'Background Epistasis':epistasis(background,mut1,mut2),
           'Positive':0,
           'Negative':0,
           'Sign':0,
           'No':0}
    
    for mut in mutations:
        if not ((mut != mut1) and (mut != mut2)):
            continue
        mutant = get_mutant(background,mut)
        ans[epistasis(mutant,mut1,mut2)] += 1

    return ans

epistasis_data = {'Background':[],
                  'Mutation 1':[],
                  'Mutation 2':[],
                  'Epistasis':[]}

for combination in nCr_index(5,2):
    
    mut1 = mutations[combination[0]]
    mut2 = mutations[combination[1]]

    for background in fitness.keys():
        epistasis_data['Background'].append(background)
        epistasis_data['Mutation 1'].append(f'{mut1[0]} at {mut1[1] + 1}')
        epistasis_data['Mutation 2'].append(f'{mut2[0]} at {mut2[1] + 1}')
        epistasis_data['Epistasis'].append(epistasis(background,mut1,mut2))

epi_df = pd.DataFrame(epistasis_data)
epi_df.to_csv('../Beta lactamase/Epistasis data.csv',index = False)

change_data = {'Background':[],
               'Mutation 1':[],
               'Mutation 2':[],
               'Background Epistasis':[],
               'Positive':[],
               'Negative':[],
               'Sign':[],
               'No':[]}

for combination in nCr_index(5,2):
    
    mut1 = mutations[combination[0]]
    mut2 = mutations[combination[1]]
    
    for background in fitness.keys():
        change_data['Background'].append(background)
        change_data['Mutation 1'].append(f'{mut1[0]} at {mut1[1] + 1}')
        change_data['Mutation 2'].append(f'{mut2[0]} at {mut2[1] + 1}')

        result = epistasis_change(background,mut1,mut2)

        for key in result.keys():
            change_data[key].append(result[key])

change_df = pd.DataFrame(change_data)
change_df.to_csv('../Beta lactamase/Epistasis change data.csv',index = False)
