from _Epistasis_functions import *

# All mutation epistasis data
path = '../Epistasis/'

mutations = generate_mutations(2)
total = len(mutations)

functional = {
    'Positive Epistasis Fraction':[],
    'Negative Epistasis Fraction':[],
    'Reciprocal Sign Epistasis Fraction':[],
    'Single Sign Epistasis Fraction':[],
    'Other Sign Epistasis Fraction':[],
    'No Epistasis Fraction':[],
    'Samples':[]
    }

non_functional = {i:[] for i in functional.keys()}
all_data = {i:[] for i in functional.keys()}

count = 0
for mutation in mutations:
    
    fun = epistasis_fraction(mutation,'functional')
    non = epistasis_fraction(mutation,'non functional')
    allp = epistasis_fraction(mutation,'all')

    for key in functional.keys():
        functional[key].append(fun[key])
        non_functional[key].append(non[key])
        all_data[key].append(allp[key])

    count += 1
    if count in [total//(100/i) for i in range(10,101,10)]:
            print(f'{round(100*count/total,0)}% mutation pairs analysed')

pd.DataFrame(functional).to_csv(f'{path}Epistasis Fractions - Functional Variants.csv',index = False)
pd.DataFrame(non_functional).to_csv(f'{path}Epistasis Fractions - Non Functional Variants.csv',index = False)
pd.DataFrame(all_data).to_csv(f'{path}Epistasis Fractions - All Variants.csv',index = False)

#Compiled epistasis data - Functional
print('\n'+'==*'*20+'\nCompiling epistasis change - functional')

save_path = '../Epistasis/Epistasis Change/'

total = len(mutations)

Positive_data = {
    
        'Background':[],
        'Background Fitness':[],
        's12 - (s1 + s2)':[],
        'Mutation 1':[],
        'Mutation 2':[],
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

total_mutations = len(generate_mutations(2))
count = 0

for mutations in generate_mutations(2):
    Pos,Neg,Rec_sign,Sin_sign,Oth_sign,No_epi = epistasis_change_data(mutations=mutations,background_nature = 'functional')
    
    name1 = f'{mutations[0][0]} at {mutations[0][1]+1}'
    name2 = f'{mutations[0][0]} at {mutations[0][1]+1}'
    
    Pos['Mutation 1'] = [name1]*len(Pos['Positive'])
    Pos['Mutation 2'] = [name2]*len(Pos['Positive'])
    
    Neg['Mutation 1'] = [name1]*len(Neg['Positive'])
    Neg['Mutation 2'] = [name2]*len(Neg['Positive'])

    Rec_sign['Mutation 1'] = [name1]*len(Rec_sign['Positive'])
    Rec_sign['Mutation 2'] = [name2]*len(Rec_sign['Positive'])

    Sin_sign['Mutation 1'] = [name1]*len(Sin_sign['Positive'])
    Sin_sign['Mutation 2'] = [name2]*len(Sin_sign['Positive'])

    Oth_sign['Mutation 1'] = [name1]*len(Oth_sign['Positive'])
    Oth_sign['Mutation 2'] = [name2]*len(Oth_sign['Positive'])

    No_epi['Mutation 1'] = [name1]*len(No_epi['Positive'])
    No_epi['Mutation 2'] = [name2]*len(No_epi['Positive'])

    for key in Pos.keys():
        Positive_data[key] += Pos[key]
        Negative_data[key] += Neg[key]
        Reciprocal_Sign_data[key] += Rec_sign[key]
        Single_Sign_data[key] += Sin_sign[key]
        Other_Sign_data[key] += Oth_sign[key]
        No_Epistasis_data[key] += No_epi[key]
    
    count += 1
    if count in [total_mutations//(100/i) for i in range(10,101,10)]:
        print(f'{round(100*count/total_mutations,0)}% mutations checked')

pos_df = pd.DataFrame(data = Positive_data)
pos_df = pos_df.sort_values(by = 'Background Fitness')
pos_df.to_csv(f'{save_path}Positive Epistasis Data - Functional.csv',index = False)

neg_df = pd.DataFrame(data = Negative_data)
neg_df = neg_df.sort_values(by = 'Background Fitness')
neg_df.to_csv(f'{save_path}Negative Epistasis Data - Functional.csv',index = False)

sin_sign_df = pd.DataFrame(data = Single_Sign_data)
sin_sign_df = sin_sign_df.sort_values(by = 'Background Fitness')
sin_sign_df.to_csv(f'{save_path}Single Sign Epistasis Data - Functional.csv',index = False)

rec_sign_df = pd.DataFrame(data = Reciprocal_Sign_data)
rec_sign_df = rec_sign_df.sort_values(by = 'Background Fitness')
rec_sign_df.to_csv(f'{save_path}Reciprocal Sign Epistasis Data - Functional.csv',index = False)

oth_sign_df = pd.DataFrame(data = Other_Sign_data)
oth_sign_df = oth_sign_df.sort_values(by = 'Background Fitness')
oth_sign_df.to_csv(f'{save_path}Other Sign Epistasis Data - Functional.csv',index = False)

no_epi_df = pd.DataFrame(data = No_Epistasis_data)
no_epi_df = no_epi_df.sort_values(by = 'Background Fitness')
no_epi_df.to_csv(f'{save_path}No Epistasis Data - Functional.csv',index = False)

#Compiled epistasis data - Non Functional
print('\n'+'==*'*20+'\nCompiling epistasis change - non functional')

total = len(mutations)

Positive_data = {
    
        'Background':[],
        'Background Fitness':[],
        's12 - (s1 + s2)':[],
        'Mutation 1':[],
        'Mutation 2':[],
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

total_mutations = len(generate_mutations(2))
count = 0

for mutations in generate_mutations(2):
    Pos,Neg,Rec_sign,Sin_sign,Oth_sign,No_epi = epistasis_change_data(mutations=mutations,background_nature = 'non functional')
    
    name1 = f'{mutations[0][0]} at {mutations[0][1]+1}'
    name2 = f'{mutations[0][0]} at {mutations[0][1]+1}'
    
    Pos['Mutation 1'] = [name1]*len(Pos['Positive'])
    Pos['Mutation 2'] = [name2]*len(Pos['Positive'])
    
    Neg['Mutation 1'] = [name1]*len(Neg['Positive'])
    Neg['Mutation 2'] = [name2]*len(Neg['Positive'])

    Rec_sign['Mutation 1'] = [name1]*len(Rec_sign['Positive'])
    Rec_sign['Mutation 2'] = [name2]*len(Rec_sign['Positive'])

    Sin_sign['Mutation 1'] = [name1]*len(Sin_sign['Positive'])
    Sin_sign['Mutation 2'] = [name2]*len(Sin_sign['Positive'])

    Oth_sign['Mutation 1'] = [name1]*len(Oth_sign['Positive'])
    Oth_sign['Mutation 2'] = [name2]*len(Oth_sign['Positive'])

    No_epi['Mutation 1'] = [name1]*len(No_epi['Positive'])
    No_epi['Mutation 2'] = [name2]*len(No_epi['Positive'])

    for key in Pos.keys():
        Positive_data[key] += Pos[key]
        Negative_data[key] += Neg[key]
        Reciprocal_Sign_data[key] += Rec_sign[key]
        Single_Sign_data[key] += Sin_sign[key]
        Other_Sign_data[key] += Oth_sign[key]
        No_Epistasis_data[key] += No_epi[key]
    
    count += 1
    if count in [total_mutations//(100/i) for i in range(10,101,10)]:
        print(f'{round(100*count/total_mutations,0)}% mutations checked')

pos_df = pd.DataFrame(data = Positive_data)
pos_df = pos_df.sort_values(by = 'Background Fitness')
pos_df.to_csv(f'{save_path}Positive Epistasis Data - Non Functional.csv',index = False)

neg_df = pd.DataFrame(data = Negative_data)
neg_df = neg_df.sort_values(by = 'Background Fitness')
neg_df.to_csv(f'{save_path}Negative Epistasis Data - Non Functional.csv',index = False)

sin_sign_df = pd.DataFrame(data = Single_Sign_data)
sin_sign_df = sin_sign_df.sort_values(by = 'Background Fitness')
sin_sign_df.to_csv(f'{save_path}Single Sign Epistasis Data - Non Functional.csv',index = False)

rec_sign_df = pd.DataFrame(data = Reciprocal_Sign_data)
rec_sign_df = rec_sign_df.sort_values(by = 'Background Fitness')
rec_sign_df.to_csv(f'{save_path}Reciprocal Sign Epistasis Data - Non Functional.csv',index = False)

oth_sign_df = pd.DataFrame(data = Other_Sign_data)
oth_sign_df = oth_sign_df.sort_values(by = 'Background Fitness')
oth_sign_df.to_csv(f'{save_path}Other Sign Epistasis Data - Non Functional.csv',index = False)

no_epi_df = pd.DataFrame(data = No_Epistasis_data)
no_epi_df = no_epi_df.sort_values(by = 'Background Fitness')
no_epi_df.to_csv(f'{save_path}No Epistasis Data - Non Functional.csv',index = False)

assert False

#Compiled epistasis data - All
print('\n'+'==*'*20+'\nCompiling epistasis change - all')

total = len(mutations)

Positive_data = {
    
        'Background':[],
        'Background Fitness':[],
        's12 - (s1 + s2)':[],
        'Mutation 1':[],
        'Mutation 2':[],
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

total_mutations = len(generate_mutations(2))
count = 0

for mutations in generate_mutations(2):
    Pos,Neg,Rec_sign,Sin_sign,Oth_sign,No_epi = epistasis_change_data(mutations=mutations,background_nature = 'all')
    
    name1 = f'{mutations[0][0]} at {mutations[0][1]+1}'
    name2 = f'{mutations[0][0]} at {mutations[0][1]+1}'
    
    Pos['Mutation 1'] = [name1]*len(Pos['Positive'])
    Pos['Mutation 2'] = [name2]*len(Pos['Positive'])
    
    Neg['Mutation 1'] = [name1]*len(Neg['Positive'])
    Neg['Mutation 2'] = [name2]*len(Neg['Positive'])

    Rec_sign['Mutation 1'] = [name1]*len(Rec_sign['Positive'])
    Rec_sign['Mutation 2'] = [name2]*len(Rec_sign['Positive'])

    Sin_sign['Mutation 1'] = [name1]*len(Sin_sign['Positive'])
    Sin_sign['Mutation 2'] = [name2]*len(Sin_sign['Positive'])

    Oth_sign['Mutation 1'] = [name1]*len(Oth_sign['Positive'])
    Oth_sign['Mutation 2'] = [name2]*len(Oth_sign['Positive'])

    No_epi['Mutation 1'] = [name1]*len(No_epi['Positive'])
    No_epi['Mutation 2'] = [name2]*len(No_epi['Positive'])

    for key in Pos.keys():
        Positive_data[key] += Pos[key]
        Negative_data[key] += Neg[key]
        Reciprocal_Sign_data[key] += Rec_sign[key]
        Single_Sign_data[key] += Sin_sign[key]
        Other_Sign_data[key] += Oth_sign[key]
        No_Epistasis_data[key] += No_epi[key]
    
    count += 1
    if count in [total_mutations//(100/i) for i in range(10,101,10)]:
        print(f'{round(100*count/total_mutations,0)}% mutations checked')

pos_df = pd.DataFrame(data = Positive_data)
pos_df = pos_df.sort_values(by = 'Background Fitness')
pos_df.to_csv(f'{save_path}Positive Epistasis Data - All.csv',index = False)

neg_df = pd.DataFrame(data = Negative_data)
neg_df = neg_df.sort_values(by = 'Background Fitness')
neg_df.to_csv(f'{save_path}Negative Epistasis Data - All.csv',index = False)

sin_sign_df = pd.DataFrame(data = Single_Sign_data)
sin_sign_df = sin_sign_df.sort_values(by = 'Background Fitness')
sin_sign_df.to_csv(f'{save_path}Single Sign Epistasis Data - All.csv',index = False)

rec_sign_df = pd.DataFrame(data = Reciprocal_Sign_data)
rec_sign_df = rec_sign_df.sort_values(by = 'Background Fitness')
rec_sign_df.to_csv(f'{save_path}Reciprocal Sign Epistasis Data - All.csv',index = False)

oth_sign_df = pd.DataFrame(data = Other_Sign_data)
oth_sign_df = oth_sign_df.sort_values(by = 'Background Fitness')
oth_sign_df.to_csv(f'{save_path}Other Sign Epistasis Data - All.csv',index = False)

no_epi_df = pd.DataFrame(data = No_Epistasis_data)
no_epi_df = no_epi_df.sort_values(by = 'Background Fitness')
no_epi_df.to_csv(f'{save_path}No Epistasis Data - All.csv',index = False)
