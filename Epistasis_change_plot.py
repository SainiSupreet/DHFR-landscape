from bio_functions import *

plt.rcParams["figure.figsize"] = (50,50)
plt.rcParams['axes.linewidth'] = 20
plt.rcParams['lines.linewidth'] = 15

plt.rcParams['axes.labelpad'] = 60

plt.rcParams['boxplot.meanprops.linewidth'] = 10
plt.rcParams['boxplot.capprops.linewidth']= 10
plt.rcParams['boxplot.meanprops.markersize'] = 10
plt.rcParams['boxplot.whiskerprops.linewidth'] = 7
plt.rcParams['boxplot.flierprops.linewidth'] = 7
plt.rcParams['boxplot.boxprops.linewidth']=7

plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.constrained_layout.h_pad']=0.5
plt.rcParams['figure.constrained_layout.w_pad']=0.5
plt.rcParams['savefig.pad_inches'] = 0.5

plt.rcParams['grid.linewidth']=3

bin_size = 0.25

min_delf = min(fitness_data.values()) - max(fitness_data.values())
max_delf = max(fitness_data.values()) - min(fitness_data.values())
range_delf = float(max_delf - min_delf)
bin_values = [min_delf+(bin_size*i) for i in range(0,int(range_delf/0.25)+1)]

epistasis_types = ['Positive','Negative','Single Sign',
                   'Reciprocal Sign','Other Sign','No Epistasis']

data_type = ['Functional','Non Functional','All']

# WRT backhround epistasis
for dtype in data_type:
    for epistasis in epistasis_types:
        ep = epistasis
        if epistasis == 'No Epistasis':
            ep = 'No'
        pos = pd.read_csv(f'../Epistasis/Epistasis Change/{ep} Epistasis Data - {dtype}.csv')
        data = []
        for i in epistasis_types:
            data.append(pos[i].values)

        bp = plt.boxplot(data,showmeans = True,meanline = True,
                         showfliers = False,
                         medianprops=dict(linestyle='-', linewidth=15,
                                          color='black'),
                         meanprops=dict(linestyle='--', linewidth=15,
                                        color='red'))
        plt.legend([bp['means'][0],bp['medians'][0]],
                   ['Mean','Median'],frameon=False,fontsize=150)
        plt.xticks(range(1,7),['PE', 'NE','SSE','RSE','OSE', 'No'],
                   fontsize = 150)
        plt.ylim(0,21)
        plt.yticks(range(0,22,3),fontsize = 165)
        plt.tick_params(axis='both', length = 45,width=10)

        plt.xlabel('Mutant epistasis',fontsize = 175)
        plt.ylabel('# Mutants',fontsize = 175)

        plt.grid()

        plt.savefig(f'../Plots/Epistasis/Epistasis Change {epistasis} Background - {dtype} Variants')
        plt.close()


# WRT mutant epistasis
for dtype in data_type:
    for i in epistasis_types:

        data = []
        for x in epistasis_types:
            ep = x
            if x == 'No Epistasis':
                ep = 'No'
            pos = pd.read_csv(f'../Epistasis/Epistasis Change/{ep} Epistasis Data - {dtype}.csv')
            data.append(pos[i].values)

        bp = plt.boxplot(data,showmeans = True,meanline = True,
                         showfliers = False,
                         medianprops=dict(linestyle='-', linewidth=15,
                                          color='black'),
                         meanprops=dict(linestyle='--', linewidth=15,
                                        color='red'))
        plt.legend([bp['means'][0],bp['medians'][0]],
                   ['Mean','Median'],frameon=False,fontsize=150)
        plt.xticks(range(1,7),['PE', 'NE','SSE','RSE','OSE', 'No'],
                   fontsize = 150)
        plt.ylim(0,21)
        plt.yticks(range(0,22,3),fontsize = 165)
        plt.tick_params(axis='both', length = 45,width=10)

        plt.xlabel('Background epistasis',fontsize = 175)
        plt.ylabel('# Mutants',fontsize = 175)

        plt.grid()

        plt.savefig(f'../Plots/Epistasis/Epistasis Change {i} mutation - {dtype} Variants')
        plt.close()

