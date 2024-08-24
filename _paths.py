from bio_functions import *

def paths_validity(variant1:str,variant2:str,
                  fitness_data = fitness_data,threshold = threshold):
    '''
    Returns a boolean result, checking if any direct path
    from variant1 to variant2 passes through only functional variants
    '''
    
    paths = generate_paths(variant1,variant2)
    
    ans = False
    for variant_list in paths:
        temp = True
        for variant in variant_list:
            if variant not in fitness_data.keys():
                continue
            elif fitness_data[variant] < threshold:
                temp = False

        if temp:
            ans = True
            break
    
    return ans

def largest_subgraph_peaks(landscape:str,variants:list,fitness_data = fitness_data,threshold = threshold):
    '''
    Checks if the peaks are connected through atleast one shortest path
    '''
    print('Started')
    reduced_variants = []
    
    for variant in variants:
        if variant not in fitness_data.keys() or fitness_data[variant] <= threshold:
            continue
        reduced_variants.append(variant)
    
    peaks = find_peaks(landscape,variants)
    print(len(peaks),'peaks found')
    
    reduced_peaks = []
    for peak in peaks:
        if fitness_data[peak] > threshold:
            reduced_peaks.append(peak)
            
    print('Peaks reduced to',len(reduced_peaks),'due to threshold')
    
    checked = [False for i in range(len(reduced_peaks))]
    
    subgraphs = []
    
    for i in range(len(reduced_peaks)):
        if checked[i]:
            continue
        
        print('New subgraph found')
        new_subgraph = []
        
        for j in range(i,len(reduced_peaks)):
            if checked[j]:
                continue
            if paths_validity(reduced_peaks[i],reduced_peaks[j]):
                new_subgraph.append(reduced_peaks[j])
                checked[j] = True
        
        print(f'Found {len(new_subgraph)} long subgraph')
        
        subgraphs.append(new_subgraph)
    
    return subgraphs

def fraction_accessibility(paths:list,fitness_data = fitness_data):
    '''
    Returns the fraction of paths accessible by darwinian evolution from a list of paths
    '''
    count = 0
    accessible = 0
    
    for path in paths:
        
        isAccessible = path_accessibility(path,fitness_data=fitness_data)
        count += 1
        
        if isAccessible:
            accessible += 1
    
    return accessible/count

def check_neighbour_accessibility(landscape_name:str,target_variant:str,distance:int,
                                  save_data = False, consider_threshold = False,path = './',
                                  fitness_data = fitness_data):
    data = {
        'Background':[],
        'Background Fitness':[],
        'Target':[],
        'Target Fitness':[],
        'Fraction Accessibility':[],
        'Number of Accessible paths':[]
    }
    
    print(f'Checking {distance} distant mutants')
    
    distant_mutants = get_n_distant_mutations('X'*9,target_variant,distance)
    
    print(f'Found {len(distant_mutants)} mutants at {distance} distance')
    
    for variant in distant_mutants:
        
        if variant not in fitness_data.keys() or (consider_threshold and fitness_data[variant] <= threshold):
            continue
        
        paths = generate_paths(variant,target_variant)
        accessibility = fraction_accessibility(paths)
        
        data['Background'].append(variant)
        data['Background Fitness'].append(fitness_data[variant])
        data['Target'].append(target_variant)
        data['Target Fitness'].append(fitness_data[target_variant])
        data['Fraction Accessibility'].append(accessibility)
        data['Number of Accessible paths'].append(int(accessibility*len(paths)))
    
    count_data = get_count_from_list(data['Number of Accessible paths'])
    
    if save_data:
        
        dataFrame = pd.DataFrame(data = data)
        dataFrame.to_csv(f'{path}{distance} distance from {target_variant}.csv',
                         index=False)
        
        count_dataFrame = pd.DataFrame(data = {'Number of Accessible paths':count_data.keys(),
                                               'Count':count_data.values()})
        count_dataFrame = count_dataFrame.sort_values(by = 'Number of Accessible paths')
        count_dataFrame.to_csv(f'{path}Count data for {distance} distance from {target_variant}.csv',
                               index=False)
        
        fig,ax = plt.subplots()
        ax.stem(count_dataFrame['Number of Accessible paths'],count_dataFrame['Count'])
        ax.set(xlabel = 'Number of Accessible paths',ylabel = 'Count',
               title = f'{distance} distance from {target_variant}')
        plt.savefig(f'{path}{distance} distance from {target_variant}')
        plt.close()
    
    return data,count_data
