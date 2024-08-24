from data import *

def get_count_from_list(data:list):
    '''
    Returns count of the number of instances a value is appearing in data
    '''
    ans = {}
    
    for val in data:
        if val in ans:
            ans[val] += 1
        else:
            ans[val] = 1
    
    return ans


def list_to_str(list1):
    '''
    Converts a list to string
    (internal function)
    '''
    ans = ''
    for value in list1:
        ans += str(value)
    
    return ans


def str_to_list(str1):
    '''
    Converts a string to list
    (internal function)
    '''
    ans = []
    for value in str(str1):
        ans.append(value)
    
    return ans


def mean(list1:list):
    '''
    Finds mean of values in list (must be int or float)
    '''
    total = 0
    for value in list1:
        total += value
    
    mean = total/len(list1)
    
    return mean


def median(list1:list):
    '''
    Finds median of values in list
    '''
    sorted_list = sorted(list1)
    
    length = len(list1)
    
    if length%2 == 0:
        ans = (sorted_list[(length//2)-1] + sorted_list[length//2])/2
    else:
        ans = sorted_list[((length+1)//2)-1]
    
    return ans


def variance(list1:list):
    '''
    Finds variance of values in list (must be int or float)
    '''
    avg = mean(list1)
    
    if len(list1) < 2:
        return 0
    
    sse = 0 # Sum of Square Error
    for value in list1:
        sse += (value-avg)**2
    
    var = sse/(len(list1)-1)
    
    return var


def find_value(list1,value):
    '''
    Returns list of indexes of value in list1
    (internal function)
    '''
    ans = []
    for i in range(len(list1)):
        if list1[i] == value:
            ans.append(i)
    return ans


def replace_instance(list1,instance,value):
    '''
    Replaces a repreating instance in a list with a
    string or list of values in the given order
    (internal function)
    '''
    instances = find_value(list1,instance)
    
    if len(instances)!=len(value):
        raise ValueError(f"{instances} Instances for {len(value)} Values")
    
    for i in range(len(instances)):
        list1[instances[i]] = value[i]
    
    return list1


def get_base_position(number,base):
    '''
    Finds the highest power of base for a integer (0<base<10)
    exmaple: highest base power for interger 16 with base 2 = 4 (2^4)
    (internal function)
    '''
    
    position = 0
    while base ** position < number:
        position += 1
        if base ** position > number:
            position -= 1
            break
    
    return position


def n_base(number,base):
    '''
    Converts an 10-base integer to any n-base (1<n<10)
    (internal function)
    '''
    
    ans = ''
    position = get_base_position(number,base)
        
    while number > 0 or position >= 0:
        new_position = get_base_position(number,base)
        
        if position == new_position:
            ans += str(number // (base**position))
            number -= int(ans[-1])*(base **position)
        else:
            ans += '0'
        position -= 1
    
    return ans


def index_to_variant(index,nucleotides = nucleotides):
    '''
    Converts 4 base index to equivalent variant variant
    (internal function)
    '''
    
    ans = ''
    
    for i in index:
        ans += nucleotides[int(i)]
    
    return ans


def remove_values(string1:str,string2:str):
    '''
    Removes instances of string 2 from string 1
    '''
    ans = ''
    
    for val in string1:
        if val not in string2:
            ans+=val
    
    return ans


def find_common_values(list1:list,list2:list):
    '''
    Returns list of common values between list1 and list2
    '''
    ans = []
    
    for val in list1:
        if val in list2:
            ans.append(val)
    
    return ans


def nCr_index(n:int,r:int):
    '''
    Finds combinations of r choose from '0…(n-1)'
    '''
    binary = [bin(i)[2:] for i in range(2**n)]
    
    filtered = []
    for val in binary:
        if list(val).count('1') == r:
            filtered.append(val)
    
    for i in range(len(filtered)):
        while len(filtered[i]) < n:
            filtered[i] = '0' + filtered[i]
    
    ans = []
    for val in filtered:
        ans.append(find_value(val,'1'))
    
    return ans


def covariance(x:list,y:list):
    '''
    Returns covariance between list of values in x and y
    '''
    if len(x) != len(y) or len(x) <= 1:
        raise ValueError('Data size invalid')
    
    x_mean = mean(x)
    y_mean = mean(y)
    
    ans = 0
    
    for i in range(len(x)):
        ans += (x[i]-x_mean)*(y[i]-y_mean)
    
    ans /= i
    
    return ans


def correlation(x:list,y:list):
    '''
    Returns correlation between values in list x and y
    '''
    covar = covariance(x,y)
    var_x = variance(x)
    var_y = variance(y)
    
    return covar/((var_x*var_y)**0.5)
