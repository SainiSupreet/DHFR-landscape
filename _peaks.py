from bio_functions import *

def get_peak_counts(n:int,save_data = False,consider_threshold = False,path = './',
                    fitness_data = fitness_data,threshold = threshold):
    '''
    Find peak counts in a 'n' point landscape
    '''
    landscapes = get_n_nucleotide_landscapes(n)
    equiv_variants = get_variants(n)
    
    number_of_landscapes = len(landscapes)
    
    print(f'Found {number_of_landscapes} landscapes and {len(equiv_variants)} variants in them')
    
    peak_count = []
    functional_variants = []
    NK_peaks = []
    
    count = 0
    for landscape in landscapes:

        functional_count = 0
        variants = [list_to_str(replace_instance(str_to_list(landscape),'X',
                                                 variant))\
                    for variant in equiv_variants]

        for variant in variants:
            if variant in fitness_data and fitness_data[variant]>threshold:
                functional_count += 1
        
        peaks = find_peaks(landscape,variants)
        
        if consider_threshold:
            
            temp_peaks = peaks.copy()
            
            peaks = []
            
            for variant in temp_peaks:
                if fitness_data[variant] > threshold:
                    peaks.append(variant)
        
        peak_count.append(len(peaks))
        functional_variants.append(functional_count)
        NK_peaks.append(round(functional_count/((3*n)+1)))
        
        count += 1
        
        if count in [number_of_landscapes//(100/i) for i in range(10,101,10)]:
            print(f'{round(100*count/number_of_landscapes,0)}% landscapes analysed')
    
    if save_data:
        
        print('Saving data')
        
        data = pd.DataFrame(
            {'Landscapes':landscapes,
             'Number of Peaks':peak_count,
             'Number of functional variants':functional_variants,
             'Predicted (NK) peaks':NK_peaks})
        
        data.to_csv(f'{path}Detailed peaks in {n} point landscapes.csv',
                    index = False)
        
        print('Data saved')
        
    return peak_count,NK_peaks

def get_peak_probability(n:int,save_data = False, consider_threshold = False,path = './',
                         fitness_data = fitness_data,threshold = threshold):
    '''
    Finds peak probability in n-point landscapes
    '''
    sequence_spaces = get_n_nucleotide_landscapes(n)
    equiv_variants = get_variants(n)
    
    number_of_variants = len(equiv_variants)
    number_of_landscapes = len(sequence_spaces)
    
    print(f'Found {number_of_landscapes} landscapes and {len(equiv_variants)} variants in them')
    
    peak_prob = []
    expected_prob = []
    
    count = 0
    
    for landscape in sequence_spaces:

        functional_count = 0
        
        variants = [list_to_str(replace_instance(str_to_list(landscape),'X',variant))\
                    for variant in equiv_variants]

        for variant in variants:
            if variant in fitness_data.keys() and fitness_data[variant]>threshold:
                functional_count += 1
        
        peaks = find_peaks(landscape,variants)
        
        if consider_threshold:
            temp_peaks = []
            
            for point in peaks:   
                if fitness_data[point] > threshold:
                    temp_peaks.append(point)
            peaks = temp_peaks
        
        peak_prob.append(len(peaks)/number_of_variants)
        expected_prob.append(round(functional_count/((3*n)+1))/(4**n))
        
        count += 1
        if count in [number_of_landscapes//(100/i) for i in range(10,101,10)]:
            print(f'{round(100*count/number_of_landscapes,0)}% landscapes analysed')
        
    if save_data:
            
        print('Saving data')
        
        data = pd.DataFrame({'Landscapes':sequence_spaces,
                             'Peak Probability':peak_prob,
                             'Expected (NK) probability':expected_prob})
        data.to_csv(f'{path}Detailed peak probability in {n} point landscapes.csv',
                    index = False)
    
        print('Data saved')
        
    return peak_prob,expected_prob
