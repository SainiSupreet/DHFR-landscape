def dfe_bins(background_fitness:list,sc:list,bin_size = 0.05,do_round = True):
    min_val = min(background_fitness)
    max_val = max(background_fitness)
    
    data = {}
    
    for i in range(len(background_fitness)):
        if do_round:
            binn = round(float(min_val + bin_size*((background_fitness[i] - min_val)//bin_size)),2)
        else:
            binn = float(min_val + bin_size*((background_fitness[i] - min_val)//bin_size))
        try: 
            try:
                data[binn].append(float(sc[i]))
            except KeyError:
                data[binn] = [float(sc[i])]
        except:
            try:
                data[binn].append(sc[i])
            except KeyError:
                data[binn] = [sc[i]]
    
    return data

def bin_data(data:list,bin_size = 0.01,max_preset = None,min_preset = None,do_round = True):
    if min_preset == None:
        min_val = min(data)
    else:
        min_val = min_preset
    
    if max_preset == None:
        max_val = max(data)
    else:
        max_val = max_preset
    
    ans = {i:0 for i in [round(float(min_val +bin_size*i),2) for i in range(0,int((max_val - min_val)//bin_size)+1,1)]}
    
    for i in range(len(data)):
        if do_round:
            binn = round(float(min_val + bin_size*int((data[i] - min_val)//bin_size)),2)
        else:
            binn = float(min_val + bin_size*int((data[i] - min_val)//bin_size))
        ans[binn] += 1
    
    return ans
