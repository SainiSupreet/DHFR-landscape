from _paths import *

target = 'AAGGAAATG'

folder = '../Peak Accessibility/'

def factorial(n):
    '''
    Returns n!
    '''
    if n > 1:
        return n*factorial(n-1)
    elif n == 1 or n == 0:
        return 1
    else:
        return 0

# Functional variants only
print(f'Checking for Functional variants')
for i in range(1,9+1):
    data = {
        'Background':[],
        'Background Fitness':[],
        'Target':[],
        'Target Fitness':[],
        'Fraction Accessibility':[],
    }
    
    print(f'Checking {i} distant mutants')
    
    distant_mutants = get_n_distant_mutations('X'*9,target,i)
    
    print(f'Found {len(distant_mutants)} mutants at {i} distance')
    
    for variant in distant_mutants:
        
        if variant not in fitness_data.keys() or fitness_data[variant] <= threshold:
            continue
            
        accessibility = fraction_accessibility(generate_paths(variant,target))
        
        data['Background'].append(variant)
        data['Background Fitness'].append(fitness_data[variant])
        data['Target'].append(target)
        data['Target Fitness'].append(fitness_data[target])
        data['Fraction Accessibility'].append(accessibility)
    
    count_data = get_count_from_list(data['Fraction Accessibility'])
    
    dataFrame = pd.DataFrame(data = data)
    dataFrame.to_csv(f'{folder}{i} distance from {target} - Functional.csv',
                     index=False)

    total = factorial(i)
    summ = sum(list(count_data.values()))
    
    count_data = pd.DataFrame(data = {'Number of Accessible Paths':[total*i for i in list(count_data.keys())],
                                      'Fraction Accessibility':count_data.keys(),
                                      'Count':count_data.values(),
                                      'PDF':[i/summ for i in list(count_data.values())]})
    count_data = count_data.sort_values(by = 'Fraction Accessibility')
    count_data.to_csv(f'{folder}Count data for {i} distance from {target} - Functional.csv',
                      index=False)

# Non-functional variants only
print(f'Checking for Non Functional variants')
for i in range(1,9+1):
    data = {
        'Background':[],
        'Background Fitness':[],
        'Target':[],
        'Target Fitness':[],
        'Fraction Accessibility':[],
    }
    
    print(f'Checking {i} distant mutants')
    
    distant_mutants = get_n_distant_mutations('X'*9,target,i)
    
    print(f'Found {len(distant_mutants)} mutants at {i} distance')
    
    for variant in distant_mutants:
        
        if variant not in fitness_data.keys() or fitness_data[variant] >= threshold:
            continue
            
        accessibility = fraction_accessibility(generate_paths(variant,target))
        
        data['Background'].append(variant)
        data['Background Fitness'].append(fitness_data[variant])
        data['Target'].append(target)
        data['Target Fitness'].append(fitness_data[target])
        data['Fraction Accessibility'].append(accessibility)
    
    count_data = get_count_from_list(data['Fraction Accessibility'])
    
    dataFrame = pd.DataFrame(data = data)
    dataFrame.to_csv(f'{folder}{i} distance from {target} - Non Functional.csv',
                     index=False)

    total = factorial(i)
    summ = sum(list(count_data.values()))
    
    count_data = pd.DataFrame(data = {'Number of Accessible Paths':[total*i for i in list(count_data.keys())],
                                      'Fraction Accessibility':count_data.keys(),
                                      'Count':count_data.values(),
                                      'PDF':[j/summ for j in list(count_data.values())]})
    count_data = count_data.sort_values(by = 'Fraction Accessibility')
    count_data.to_csv(f'{folder}Count data for {i} distance from {target} - Non Functional.csv',
                      index=False)

# All variants
print(f'Checking for All variants')
for i in range(1,9+1):
    data = {
        'Background':[],
        'Background Fitness':[],
        'Target':[],
        'Target Fitness':[],
        'Fraction Accessibility':[],
    }
    
    print(f'Checking {i} distant mutants')
    
    distant_mutants = get_n_distant_mutations('X'*9,target,i)
    
    print(f'Found {len(distant_mutants)} mutants at {i} distance')
    
    for variant in distant_mutants:
        
        if variant not in fitness_data.keys():
            continue
            
        accessibility = fraction_accessibility(generate_paths(variant,target))
        
        data['Background'].append(variant)
        data['Background Fitness'].append(fitness_data[variant])
        data['Target'].append(target)
        data['Target Fitness'].append(fitness_data[target])
        data['Fraction Accessibility'].append(accessibility)
    
    count_data = get_count_from_list(data['Fraction Accessibility'])
    
    dataFrame = pd.DataFrame(data = data)
    dataFrame.to_csv(f'{folder}{i} distance from {target}.csv',
                     index=False)

    total = factorial(i)
    summ = sum(list(count_data.values()))
    
    count_data = pd.DataFrame(data = {'Number of Accessible Paths':[total*i for i in list(count_data.keys())],
                                      'Fraction Accessibility':count_data.keys(),
                                      'Count':count_data.values(),
                                      'PDF':[i/summ for i in list(count_data.values())]})
    count_data = count_data.sort_values(by = 'Fraction Accessibility')
    count_data.to_csv(f'{folder}Count data for {i} distance from {target} - All.csv',
                      index=False)
