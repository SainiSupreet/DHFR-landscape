from math_functions import *

def H_distance(variant1,variant2):
    '''
    Finds Hamming Distance between two variants
    '''
    
    if len(variant1) != len(variant2): #If variants are of different lengths
        raise ValueError(f"Variants {variant1} and {variant2} are of different length")
    
    ans = 0
    for i in range(len(variant1)):
        if variant1[i] != variant2[i]:
            ans += 1
    
    return ans


def get_variants(order:int):
    '''
    Generate all possible variants of a given order (number of involved nucleotides)
    '''
    
    variants = []
    
    # If number of nucleotides is < 1, return the empty list
    if order < 1:
        return variants
    
    # Considering all possible 4^n combinations
    for i in range(4**order):
        
        # Creating a 4 base integer to represent the genotype
        # (0,1,2,3 represent A,T,G,C respectively)
        variant_index = n_base(i,4)
        
        # Ensuring the resulting genotype is of 'n' length
        while len(variant_index) < order:
            variant_index = '0'+variant_index
        
        variants.append(index_to_variant(variant_index))
    
    return variants


def get_floating_positions(order:int):
    '''
    Generate all floating positions of the given order
    (internal function)
    '''
    
    floating_positions = []
    
    # Considering all possible 2^9 combinations
    for i in range(1,2**9):
        
        # Generating a binary number
        # ('0' corresporing to fixed positions and '1' corresponding to floating position)
        floating_index = bin(i)[2:]
        
        # Seperating all landscapes with desierd number of '1'
        if floating_index.count('1')!= order:
            continue
        
        # Ensuring all results are of length 9
        while len(floating_index) < 9:
            floating_index = '0' + floating_index
        
        floating_positions.append(floating_index.replace('1','X'))
    
    return floating_positions


def get_n_nucleotide_landscapes(n:int):
    '''
    Generates all possible landscapes of 'n' nucleotides
    '''
    
    floating_positions = get_floating_positions(n) # All permutations of floating positions
    fixed_values = get_variants(9-n) # All permutations of fixed points
    
    ans = []
    
    # For a given orientation of floating positions, fill in the fixed positions
    for floating in floating_positions:
        
        if len(fixed_values) > 0:
            for filler in fixed_values:
                
                # Filling in the fixed position
                floating_list = replace_instance(str_to_list(floating),'0',filler)
                ans.append(list_to_str(floating_list))
                
        else: # In case no fixed position exist, return the floating position list (case for n = 9)
            ans.append(floating)
    
    return ans


def check_availability(landscape_name:str,variant:str):
    '''
    Checks if the variant is in the landscape
    '''
    
    if len(landscape_name) != len(variant):
        return False
    
    loci = find_value(landscape_name,'X')
    
    ans = True
    
    for i in range(len(landscape_name)):
        if i in loci:
            continue
        else:
            if landscape_name[i] != variant[i]:
                ans = False
                break
    
    return ans


def open_landscape(landscape:str):
    '''
    Returns all variants in the landscape
    '''
    
    order = landscape.count('X') # Finding the number of nucleotides contributing to the landscape
    variants = get_variants(order) # Find all permutations of points in the landscape of equivalent nucleotides
    
    # Generating list of all variants in the landscape by replacing the floating nucleotides with 'variants' list
    ans = [list_to_str(replace_instance(str_to_list(landscape),'X',variant)) 
           for variant in variants]
    
    return ans


def find_neighbours(landscape:str, variant:str,nucleotides = nucleotides):
    '''
    Finds neighbours of the variant in a landscape
    '''
    
    floating_locations = find_value(landscape,'X')
        
    ans = []
    
    for index in floating_locations:
        for nucleotide in nucleotides:
            if variant[index] != nucleotide:
                
                variant_list = str_to_list(variant)
                variant_list[index] = nucleotide
                
                ans.append(list_to_str(variant_list))
                
                del variant_list
            
    
    return ans


def get_n_distant_mutations(landscape_name:str,variant:str,n:int):
    '''
    Returns list of mutations at n distance from the variant in the landscape
    '''
    if n == 0:
        return [variant]
    
    if landscape_name.count('X') < n or n < 0:
        raise ValueError(f'Cannot find {n} distant neighbours in the landscape')
    elif not check_availability(landscape_name,variant):
        raise ValueError('Variant not in the landscape')
    
    loci = find_value(landscape_name,'X')
    
    n_mutants = []
    
    # list of sites on which two mutations can act on
    sites = []
    for site in nCr_index(len(loci),n):
            sites.append([loci[int(i)] for i in site])
    
    for site in sites:
        for combination in get_variants(n):
            
            if True in [variant[site[i]] == combination[i] for i in range(n)]:
                continue
            
            n_mutant = str_to_list(variant)
            for i in range(n):
                n_mutant[site[i]] = combination[i]
            
            n_mutants.append(list_to_str(n_mutant))
    
    return n_mutants


def involvedVariants(variant1:str,variant2:str):
    '''
    Returns a complex list of all variants involed in direct paths from variant 1 to variant 2
    (internal function)
    '''
    
    if len(variant1) != len(variant2):
        raise ValueError('Sequences are not of equal length')
        
    distance = H_distance(variant1,variant2)
    
    landscape_name = ''
    for i in range(len(variant1)):
        if variant1[i] == variant2[i]:
            landscape_name += variant1[i]
        else:
            landscape_name += 'X'
    
    ans = []
    
    for h in range(distance+1):
        
        if h < distance/2:
            
            list1 = get_n_distant_mutations(landscape_name,variant1,h)
            sub_ans = []
            
            for val in list1:
                if H_distance(val,variant2) == distance - h:
                    sub_ans.append(val)
            
            ans.append(sub_ans)
        
        else:
            
            list1 = get_n_distant_mutations(landscape_name,variant2,distance-h)
            sub_ans = []
            
            for val in list1:
                if H_distance(val,variant1) == h:
                    sub_ans.append(val)
            
            ans.append(sub_ans)
    
    return ans


def make_path_list(variants:list,explored_paths:list = []):
    '''
    Generates list of all paths using the 'involved variants' list
    (internal function)
    '''
    
    ans = explored_paths.copy()
    
    if len(ans) == 0: #First step
        ans.append(variants[0])
    
    
    step = len(ans[0]) #Number of steps already calculated
    
    if step == len(variants): #If all paths are fully explored, return final answer
        return ans
    
    new_ans = [] #Since lists are mutable, new list is created to avoid mutability errors
    
    for path in ans:
        
        # List of possible next steps in interated 'path'
        next_steps = []
        
        for variant in variants[step]:
            if H_distance(variant,path[-1]) == 1:
                next_steps.append(variant)
        
        # Creating branches of each path
        for variant in next_steps:
            
            path_copy = path.copy()
            path_copy.append(variant)
            
            new_ans.append(path_copy)
    
    # Reccursively solve to find all possible paths
    return make_path_list(variants,new_ans)

def generate_paths(variant1:str,variant2:str):
    '''
    Generates list of direct paths from a variant1 to variant2 in the landscape
    '''
    
    involved_variants = involvedVariants(variant1,variant2) # Finding involved variants
    return make_path_list(involved_variants) # Generating paths from the involved variants

def find_peaks(landscape:str,landscape_list:list):
    '''
    Returns list of peaks in the given landscape
    '''
    
    peaks = [] # list of peaks
    
    for point in landscape_list:
        
        # If data for point is not availaible, it is not a peak
        if point not in fitness_data.keys():
            continue
        
        neighbours = find_neighbours(landscape,point) # list of neighbours of point in the landscape
        
        isPeak = True
        
        for neighbour in neighbours:
            
            # If data for point is not availaible, it cannot be compared
            if neighbour not in fitness_data.keys():
                continue
            
            # If neighbour is a peak or has a fitness greater than the current point,
            # the point can not be a peak
            if (neighbour in peaks) or (fitness_data[neighbour] > fitness_data[point]):
                isPeak = False
                break
        
        if isPeak:
            peaks.append(point)
        
    return peaks


def main_subgraph(landscape:str,landscape_list:list,threshold = threshold):
    '''
    List of functional variants
    '''
    peaks = find_peaks(landscape,landscape_list)
    max_peak = peaks[0]
    
    for peak in peaks:
        if fitness_data[peak] > fitness_data[max_peak]:
            max_peak = peak
    
    large_subgraph = [max_peak]
    variants = open_landscape(landscape)
    
    wave = 0
    
    change = True
    while change == True:
        
        wave += 1
        count = 0
        
        change = False
        for variant in variants:
            if variant in large_subgraph:
                continue
                
            for neighbour in find_neighbours(landscape,variant):
                if neighbour not in fitness_data.keys():
                    continue
                if variant not in fitness_data.keys() and fitness_data[neighbour] < threshold:
                    continue
                if (fitness_data[neighbour] > threshold or fitness_data[variant] > threshold)\
                and neighbour in large_subgraph:
                    change = True
                    large_subgraph.append(variant)
                    count += 1
                    break
                    
        print(f'Wave {wave}, Newly found points = {count}')
    
    return large_subgraph

def find_common_neighbours(landscape_name:str,variant1:str,variant2:str):
    '''
    Returns list of common neighbours between variant1 and variant2 in the landscape
    '''
    neighbours1 = find_neighbours(landscape=landscape_name,variant=variant1)
    neighbours2 = find_neighbours(landscape=landscape_name,variant=variant2)
    
    ans = []
    
    for val in neighbours1:
        if val in neighbours2:
            ans.append(val)
    
    return ans


def path_accessibility(path:list,fitness_data = fitness_data):
    '''
    Returns weather a path is accessible by adaptive mutations
    '''
    isAccessible = True
    
    fitness = -100
    
    for variant in path:
        if (variant not in fitness_data.keys()) or (fitness_data[variant] <= fitness):
            isAccessible = False
            break
        
        fitness = fitness_data[variant]
    
    return isAccessible


def generate_mutations(distance = 1,base_pairs = 9,nucleotides = nucleotides):
    '''
    Generates all mutation combinations in (base,loci) format
    '''
    if distance > base_pairs:
        raise ValueError(f'Cannot find {distance} distant mutations in a {base_pairs} long variant')
    
    mutations_loci = nCr_index(base_pairs,distance)
    variants = get_variants(distance)
    
    mutations = []
    
    for loci in mutations_loci:
        for bases in variants:
            mutation = []
            for i in range(distance):
                mutation.append((bases[i],loci[i]))
            
            mutations.append(mutation)
    
    return mutations


def mutation_feasibility(variant:str,mutation:tuple):
    '''
    Checks if mutation is feasible on the given variant
    '''
    try:
        return mutation[0] != variant[mutation[1]]
    except IndexError:
        return False


def get_mutant(variant:str,mutation:tuple):
    '''
    Get mutant sequence of the given variant
    '''
    if len(variant) < mutation[1] + 1:
        raise ValueError(f'Variant cannot have a mutation on locus {mutation[1] + 1}')
    list_equiv = list(variant)
    list_equiv[mutation[1]] = mutation[0]
    
    return list_to_str(list_equiv)


def check_epistasis_nature(background:str,mutations:list,fitness_data = fitness_data,
                           no_epi_margin = no_epi_margin):
    '''
    Checks if the mutation pairs result in Positive, Negative or Sign Epistasis on the given background sequence
    '''
    
    if len(mutations) != 2 or mutations[0][1] == mutations[1][1]\
    or not mutation_feasibility(background,mutations[0]) or not mutation_feasibility(background,mutations[1]):
        raise ValueError('Mutation Pairs are invalid')
    
    mutant1 = get_mutant(background,mutations[0])
    mutant2 = get_mutant(background,mutations[1])
    double_mutant = get_mutant(mutant1,mutations[1])
    
    if background not in fitness_data.keys()\
    or mutant1 not in fitness_data.keys()\
    or mutant2 not in fitness_data.keys()\
    or double_mutant not in fitness_data.keys():
        return {'Nature':'Unknown'}
    
    s1 = fitness_data[mutant1] - fitness_data[background]
    s2 = fitness_data[mutant2] - fitness_data[background]
    s12 = fitness_data[double_mutant] - fitness_data[background]
    
    net_sc = s1 + s2
    
    if abs(net_sc - s12) < no_epi_margin:
        return {'Nature':'No Epistasis','s12 - (s1 + s2)':s12 - net_sc}
    elif net_sc*s12 < 0:
        if s1 < 0 and s2 < 0 and s12 > 0:
            return {'Nature':'Reciprocal Sign','s12 - (s1 + s2)':s12 - net_sc}
        elif (s1 < 0 or s2 < 0) and s12 > 0:
            return {'Nature':'Single Sign','s12 - (s1 + s2)':s12 - net_sc}
        else:
            return {'Nature':'Other Sign','s12 - (s1 + s2)':s12 - net_sc}
    elif net_sc < s12:
        return {'Nature':'Positive','s12 - (s1 + s2)':s12 - net_sc}
    elif net_sc > s12:
        return {'Nature':'Negative','s12 - (s1 + s2)':s12 - net_sc}
