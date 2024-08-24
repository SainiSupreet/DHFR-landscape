from bio_functions import *

mutations = []

for mut in generate_mutations(1):
    mutations.append(mut[0])

types = ['functional','non functional','all']
neutral_threshold = 0.05

for typ in types:
    
    for mut in mutations:

        background_pairs = []
        for var in fitness_data.keys():

            if not mutation_feasibility(var,mut):
                continue
            
            fitness = fitness_data[var]

            if (typ == 'functional' and fitness < threshold) or\
               (typ == 'non functional' and fitness > threshold):
                continue
            
            mutant = get_mutant(var,mut)

            try:
                mut_fitness = fitness_data[mutant]
            except KeyError:
                continue

            if abs(mut_fitness - fitness) > neutral_threshold:
                continue

            background_pairs.append((var,mutant))

        data = {'Background':[],
                'Background Fitness':[],
                'Mutant':[],
                'Mutant Fitness':[],
                'New Mutation':[],
                'SC 1':[],
                'SC 2':[],
                'delta SC':[],
                'relative delta SC':[]}

        for bp in background_pairs:
            bf = fitness_data[bp[0]]
            mf = fitness_data[bp[1]]

            for mut2 in mutations:
                
                if not (mutation_feasibility(bp[0],mut2) and mutation_feasibility(bp[1],mut2)):
                    continue

                m1 = get_mutant(bp[0],mut2)
                m2 = get_mutant(bp[1],mut2)

                try:
                    f1 = fitness_data[m1]
                    f2 = fitness_data[m2]
                except KeyError:
                    continue
                
                data['Background'].append(bp[0])
                data['Background Fitness'].append(bf)
                data['Mutant'].append(bp[1])
                data['Mutant Fitness'].append(mf)
                data['New Mutation'].append(f'{mut2[0]} at {mut2[1]+1}')

                sc1 = f1 - bf
                sc2 = f2 - mf
                del_sc = sc2 - sc1
                
                data['SC 1'].append(sc1)
                data['SC 2'].append(sc2)

                data['delta SC'].append(del_sc)
                data['relative delta SC'].append(del_sc/abs(sc1))

        df = pd.DataFrame(data)
        df = df.sort_values(by = 'Background Fitness')
        df.to_csv(f'../Mutation effect on evolvability/{typ} variants/{mut[0]} at {mut[1]+1}.csv',
                  index = False)
