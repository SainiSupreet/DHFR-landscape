from bio_functions import *

def find_bin(bins:list,value:float):
    '''
    Internal function
    '''
    
    for i in range(len(bins)):
        if bins[i] > value:
            i -= 1
            break

    if i == -1:
        raise ValueError(f'bins {bins} too large for given value {value}')

    return i

def epistasis_fraction(mutations,background_nature = ['functional','non functional','all'],
                       fitness_data=fitness_data,threshold = threshold,
                       no_epi_margin = no_epi_margin):
    '''
    Returns Compiled data of the fraction of Benefecial, Deleterious and 
    Sign Epistasis present in the mutation pairs
    '''
    
    if len(mutations) != 2:
        raise ValueError('Exactly 2 mutations are required')

    if background_nature not in ['functional','non functional','all']:
        raise ValueError(f'Background nature {background_nature} not valid')
    
    Backgrounds_checked = 0
    Positive_count = 0
    Negative_count = 0
    Reciprocal_Sign_count = 0
    Single_Sign_count = 0
    Other_Sign_count = 0
    No_epi_count = 0
    
    for background in fitness_data.keys():

        background_fitness = fitness_data[background]
        if (background_fitness < threshold and background_nature == 'functional') or\
           (background_fitness > threshold and background_nature == 'non functional'):
            continue

        # Skip if background has same base on mutation locus
        if background[mutations[0][1]] == mutations[0][0] or\
           background[mutations[1][1]] == mutations[1][0]:
            continue
        
        try:
            mutant1_fitness = fitness_data[get_mutant(background,mutations[0])]
            mutant2_fitness = fitness_data[get_mutant(background,mutations[1])]
            double_mutant_fitness = fitness_data[get_mutant(get_mutant(background,mutations[0]),
                                                            mutations[1])]
            
        except KeyError:
            # Skip if fitness of mutant is not known
            continue
        
        epi_data = check_epistasis_nature(background,mutations)
        
        if epi_data['Nature'] == 'No Epistasis':
            No_epi_count += 1
        elif epi_data['Nature'] == 'Reciprocal Sign':
            Reciprocal_Sign_count += 1
        elif epi_data['Nature'] == 'Single Sign':
            Single_Sign_count += 1
        elif epi_data['Nature'] == 'Other Sign':
            Other_Sign_count += 1
        elif epi_data['Nature'] == 'Positive':
            Positive_count += 1
        elif epi_data['Nature'] == 'Negative':
            Negative_count += 1
            
        Backgrounds_checked += 1
    
    ans = {
        'Positive Epistasis Fraction':Positive_count/Backgrounds_checked,
        'Negative Epistasis Fraction':Negative_count/Backgrounds_checked,
        'Reciprocal Sign Epistasis Fraction':Reciprocal_Sign_count/Backgrounds_checked,
        'Single Sign Epistasis Fraction':Single_Sign_count/Backgrounds_checked,
        'Other Sign Epistasis Fraction':Other_Sign_count/Backgrounds_checked,
        'No Epistasis Fraction':No_epi_count/Backgrounds_checked,
        'Samples':Backgrounds_checked
          }
    return ans

def get_epistasis_data(mutations,background_nature = ['functional','non functional','all'],
                       save_data = False,path = './',
                       fitness_data = fitness_data,threshold = threshold,
                       no_epi_margin = no_epi_margin):
    '''
    Returns detailed data of epistasis present in all backgrounds for a given mutation pair
    '''
    if len(mutations) != 2:
        raise ValueError('Exactly 2 mutations are required')

    if background_nature not in ['functional','non functional','all']:
        raise ValueError(f'Background nature {background_nature} not valid')

    if save_data:
        file_name = ''
        for mutation in mutations:
            if file_name != '':
                           file_name += ' & '
                
            file_name += f'{mutation[0]} at {mutation[1]+1}'
    
    data = {
        'Background':[],
        'Background Fitness':[],
        'Mutant 1':[],
        'Mutant 1 Fitness':[],
        's1':[],
        'Mutant 2':[],
        'Mutant 2 Fitness':[],
        's2':[],
        'Double Mutant':[],
        'Double Mutant Fitness':[],
        's12':[],
        's12 - (s1 + s2)':[],
        'Epistasis':[]
    }
    
    for background in fitness_data.keys():

        background_fitness = fitness_data[background]
        if (background_fitness < threshold and background_nature == 'functional') or\
           (background_fitness > threshold and background_nature == 'non functional'):
            continue
        
        if background[mutations[0][1]] == mutations[0][0] or\
           background[mutations[1][1]] == mutations[1][0]:
            continue
        
        try:
            mutant1_fitness = fitness_data[get_mutant(background,mutations[0])]
            mutant2_fitness = fitness_data[get_mutant(background,mutations[1])]
            double_mutant_fitness = fitness_data[get_mutant(get_mutant(background,mutations[0]),
                                                            mutations[1])]
            
        except KeyError:
            # Skip if fitness of mutant is not known
            continue
        
        mutants = [get_mutant(background,i) for i in mutations]
        double_mutant = get_mutant(get_mutant(background,mutations[0]),mutations[1])
        
        s1 = mutant1_fitness - background_fitness
        s2 = mutant2_fitness - background_fitness
        s12 = double_mutant_fitness - background_fitness
        
        net_sc = s1+s2
        
        data['Background'].append(background)
        data['Background Fitness'].append(background_fitness)
        
        data['Mutant 1'].append(mutants[0])
        data['Mutant 1 Fitness'].append(mutant1_fitness)
        data['s1'].append(s1)
        
        data['Mutant 2'].append(mutants[1])
        data['Mutant 2 Fitness'].append(mutant2_fitness)
        data['s2'].append(s2)
        
        data['Double Mutant'].append(double_mutant)
        data['Double Mutant Fitness'].append(double_mutant_fitness)
        data['s12'].append(s12)
        
        epi_data = check_epistasis_nature(background,mutations)
        data['s12 - (s1 + s2)'].append(epi_data['s12 - (s1 + s2)'])
        data['Epistasis'].append(epi_data['Nature'])
        
    if save_data:
        dataFrame = pd.DataFrame(data = data)
        dataFrame.to_csv(f'{path+file_name}.csv',index = False)
    
    return data

def epistasis_change_due_to_mutation(background:str,mutations:list,
                                     del_f_bins = []):
    '''
    Returns the nature of epistasis in the neighbours of the given background
    '''
    indices = list(range(len(background)))
    indices.remove(mutations[0][1])
    indices.remove(mutations[1][1])

    bin_data = True
    if len(del_f_bins)==0:
        bin_data = False
    
    equiv_seq_space = list_to_str(['X' if i in indices else background[i] for i in range(len(background))])
    
    neighbours = find_neighbours(equiv_seq_space,background)
    
    data = {
        'Positive':0,
        'Negative':0,
        'Reciprocal Sign':0,
        'Single Sign':0,
        'Other Sign':0,
        'No Epistasis':0,
        'Unknown':0
    }

    fitness_bins = {}
    
    for variant in neighbours:

        if bin_data:
            try:
                binn = find_bin(del_f_bins,fitness_data[variant]-fitness_data[background])
            except KeyError:
                continue

            nature = check_epistasis_nature(variant,mutations)
            try:
                fitness_bins[del_f_bins[binn]][nature['Nature']] += 1
            except KeyError:
                fitness_bins[del_f_bins[binn]] = {j:0 for j in data.keys()}
                fitness_bins[del_f_bins[binn]][nature['Nature']] += 1

        else:
            nature = check_epistasis_nature(variant,mutations)
            data[nature['Nature']] += 1

    if bin_data:
        return fitness_bins
    
    return data

def epistasis_change_data(mutations:tuple,background_nature = ['functional','non functional','all']):
    '''
    Returns Epistasis change data
    '''
    
    if len(mutations) != 2:
        raise ValueError('Exactly 2 mutations are required')

    if background_nature not in ['functional','non functional','all']:
        raise ValueError(f'Background nature {background_nature} not valid')
    
    Positive_data = {
        'Background':[],
        'Background Fitness':[],
        's12 - (s1 + s2)':[],
        'Positive':[],
        'Reciprocal Sign':[],
        'Single Sign':[],
        'Other Sign':[],
        'Negative':[],
        'No Epistasis':[],
        'Unknown':[]
    }
    
    Negative_data = {i:[] for i in Positive_data.keys()}
    Reciprocal_Sign_data = {i:[] for i in Positive_data.keys()}
    Single_Sign_data = {i:[] for i in Positive_data.keys()}
    Other_Sign_data = {i:[] for i in Positive_data.keys()}
    No_Epistasis_data = {i:[] for i in Positive_data.keys()}
    
    for background in fitness_data.keys():
        
        background_fitness = fitness_data[background]
        if (background_fitness < threshold and background_nature == 'functional') or\
           (background_fitness > threshold and background_nature == 'non functional'):
            continue
        
        if background[mutations[0][1]] == mutations[0][0] or\
           background[mutations[1][1]] == mutations[1][0]:
            continue
        
        try:
            mutant1_fitness = fitness_data[get_mutant(background,mutations[0])]
            mutant2_fitness = fitness_data[get_mutant(background,mutations[1])]
            double_mutant_fitness = fitness_data[get_mutant(get_mutant(background,mutations[0]),
                                                            mutations[1])]
            
        except KeyError:
            # Skip if fitness of mutant is not known
            continue
        
        epistasis_nature = check_epistasis_nature(background,mutations)
        
        if epistasis_nature['Nature'] == 'Positive':
            
            change_data = epistasis_change_due_to_mutation(background,mutations)
            
            Positive_data['Background'].append(background)
            Positive_data['Background Fitness'].append(background_fitness)
            Positive_data['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])

            for key in change_data.keys():
                Positive_data[key].append(change_data[key])
        
        elif epistasis_nature['Nature'] == 'Negative':
            
            change_data = epistasis_change_due_to_mutation(background,mutations)
            
            Negative_data['Background'].append(background)
            Negative_data['Background Fitness'].append(background_fitness)
            Negative_data['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
            
            for key in change_data.keys():
                Negative_data[key].append(change_data[key])
        
        elif epistasis_nature['Nature'] == 'Reciprocal Sign':
            
            change_data = epistasis_change_due_to_mutation(background,mutations)
            
            Reciprocal_Sign_data['Background'].append(background)
            Reciprocal_Sign_data['Background Fitness'].append(background_fitness)
            Reciprocal_Sign_data['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
            
            for key in change_data.keys():
                Reciprocal_Sign_data[key].append(change_data[key])
        
        elif epistasis_nature['Nature'] == 'Single Sign':
            
            change_data = epistasis_change_due_to_mutation(background,mutations)
            
            Single_Sign_data['Background'].append(background)
            Single_Sign_data['Background Fitness'].append(background_fitness)
            Single_Sign_data['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
            
            for key in change_data.keys():
                Single_Sign_data[key].append(change_data[key])
        
        elif epistasis_nature['Nature'] == 'Other Sign':
            
            change_data = epistasis_change_due_to_mutation(background,mutations)
            
            Other_Sign_data['Background'].append(background)
            Other_Sign_data['Background Fitness'].append(background_fitness)
            Other_Sign_data['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
            
            for key in change_data.keys():
                Other_Sign_data[key].append(change_data[key])
        
        elif epistasis_nature['Nature'] == 'No Epistasis':
            
            change_data = epistasis_change_due_to_mutation(background,mutations)
            
            No_Epistasis_data['Background'].append(background)
            No_Epistasis_data['Background Fitness'].append(background_fitness)
            No_Epistasis_data['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
            
            for key in change_data.keys():
                No_Epistasis_data[key].append(change_data[key])
        
        else:
            raise ValueError(f'No test condition for {epistasis_nature["Nature"]} found')
    
    return Positive_data,Negative_data,Reciprocal_Sign_data,Single_Sign_data,Other_Sign_data,No_Epistasis_data

def epistasis_change_del_f(mutations:tuple,
                           background_nature = ['functional',
                                                'non functional',
                                                'all'],
                           del_f_bins = []):
    '''
    Returns Epistasis change data
    '''
    
    if len(mutations) != 2:
        raise ValueError('Exactly 2 mutations are required')

    if background_nature not in ['functional','non functional','all']:
        raise ValueError(f'Background nature {background_nature} not valid')

    if len(del_f_bins) == 0:
        raise ValueError('Fitness bins must not be empty')

    ref_dict = {
        'Background':[],
        'Background Fitness':[],
        's12 - (s1 + s2)':[],
        'Positive':[],
        'Reciprocal Sign':[],
        'Single Sign':[],
        'Other Sign':[],
        'Negative':[],
        'No Epistasis':[],
        'Unknown':[]
    }
    
    Positive_data = {i:{j:[] for j in ref_dict.keys()} for i in del_f_bins}
    Negative_data = {i:{j:[] for j in ref_dict.keys()} for i in del_f_bins}
    Reciprocal_Sign_data = {i:{j:[] for j in ref_dict.keys()} for i in del_f_bins}
    Single_Sign_data = {i:{j:[] for j in ref_dict.keys()} for i in del_f_bins}
    Other_Sign_data = {i:{j:[] for j in ref_dict.keys()} for i in del_f_bins}
    No_Epistasis_data = {i:{j:[] for j in ref_dict.keys()} for i in del_f_bins}
    
    for background in fitness_data.keys():
        
        background_fitness = fitness_data[background]
        if (background_fitness < threshold and background_nature == 'functional') or\
           (background_fitness > threshold and background_nature == 'non functional'):
            continue
        
        if background[mutations[0][1]] == mutations[0][0] or\
           background[mutations[1][1]] == mutations[1][0]:
            continue
        
        try:
            mutant1_fitness = fitness_data[get_mutant(background,mutations[0])]
            mutant2_fitness = fitness_data[get_mutant(background,mutations[1])]
            double_mutant_fitness = fitness_data[get_mutant(get_mutant(background,mutations[0]),
                                                            mutations[1])]
            
        except KeyError:
            # Skip if fitness of mutant is not known
            continue
        
        epistasis_nature = check_epistasis_nature(background,mutations)
        change_data = epistasis_change_due_to_mutation(background,mutations,
                                                        del_f_bins)
        
        if epistasis_nature['Nature'] == 'Positive':
            
            for binn in change_data.keys():
                Positive_data[binn]['Background'].append(background)
                Positive_data[binn]['Background Fitness'].append(background_fitness)
                Positive_data[binn]['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
                for nature in change_data[binn].keys():
                    Positive_data[binn][nature].append(change_data[binn][nature])
        
        elif epistasis_nature['Nature'] == 'Negative':

            for binn in change_data.keys():
                Negative_data[binn]['Background'].append(background)
                Negative_data[binn]['Background Fitness'].append(background_fitness)
                Negative_data[binn]['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
                for nature in change_data[binn].keys():
                    Negative_data[binn][nature].append(change_data[binn][nature])
        
        elif epistasis_nature['Nature'] == 'Reciprocal Sign':
            
            for binn in change_data.keys():
                Reciprocal_Sign_data[binn]['Background'].append(background)
                Reciprocal_Sign_data[binn]['Background Fitness'].append(background_fitness)
                Reciprocal_Sign_data[binn]['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
                for nature in change_data[binn].keys():
                    Reciprocal_Sign_data[binn][nature].append(change_data[binn][nature])
        
        elif epistasis_nature['Nature'] == 'Single Sign':
            
            for binn in change_data.keys():
                Single_Sign_data[binn]['Background'].append(background)
                Single_Sign_data[binn]['Background Fitness'].append(background_fitness)
                Single_Sign_data[binn]['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
                for nature in change_data[binn].keys():
                    Single_Sign_data[binn][nature].append(change_data[binn][nature])
        
        elif epistasis_nature['Nature'] == 'Other Sign':
            
            for binn in change_data.keys():
                Other_Sign_data[binn]['Background'].append(background)
                Other_Sign_data[binn]['Background Fitness'].append(background_fitness)
                Other_Sign_data[binn]['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
                for nature in change_data[binn].keys():
                    Other_Sign_data[binn][nature].append(change_data[binn][nature])
        
        elif epistasis_nature['Nature'] == 'No Epistasis':
            
            for binn in change_data.keys():
                No_Epistasis_data[binn]['Background'].append(background)
                No_Epistasis_data[binn]['Background Fitness'].append(background_fitness)
                No_Epistasis_data[binn]['s12 - (s1 + s2)'].append(epistasis_nature['s12 - (s1 + s2)'])
                for nature in change_data[binn].keys():
                    No_Epistasis_data[binn][nature].append(change_data[binn][nature])
        
        else:
            raise ValueError(f'No test condition for {epistasis_nature["Nature"]} found')
        
    return Positive_data,Negative_data,Reciprocal_Sign_data,Single_Sign_data,Other_Sign_data,No_Epistasis_data
