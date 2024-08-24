from bio_functions import *

plt.rcParams["figure.figsize"] = (100,50)
plt.rcParams['axes.linewidth'] = 20
plt.rcParams['lines.linewidth'] = 15

plt.rcParams['axes.labelpad'] = 30

plt.rcParams['boxplot.meanprops.linewidth'] = 10
plt.rcParams['boxplot.capprops.linewidth']= 10
plt.rcParams['boxplot.meanprops.markersize'] = 10
plt.rcParams['boxplot.whiskerprops.linewidth'] = 7
plt.rcParams['boxplot.flierprops.linewidth'] = 7
plt.rcParams['boxplot.boxprops.linewidth']=7

plt.rcParams['xtick.major.size']=50
plt.rcParams['xtick.major.width']=10
plt.rcParams['ytick.major.size']=50
plt.rcParams['ytick.major.width']=10

plt.rcParams['grid.linewidth']=3

plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.constrained_layout.h_pad']=0.5
plt.rcParams['figure.constrained_layout.w_pad']=0.5
plt.rcParams['savefig.pad_inches'] = 0.5

types = ['functional','non functional','all']
mutations = [i[0] for i in generate_mutations(1)]

for typ in types:
    data = []
    for mut in mutations:
        file_name = f'{mut[0]} at {mut[1]+1}'
        file = pd.read_csv(f'../Mutation effect on evolvability/{typ} variants/{file_name}.csv')
        data.append(file['relative delta SC'].values)

    bp = plt.boxplot(data, showfliers=False,
                medianprops=dict(linestyle='-', linewidth=15, color='red'))

    p = plt.plot((0.5,36.5),(0,0),linestyle = 'dashed',color = 'black',alpha = 0.4)

    plt.xlabel('Mutation',fontsize = 175)
    plt.ylabel('Relative increase in fitness effect',fontsize = 175)
    plt.xticks(range(1,37),[f'{i[0]}\n{i[1]+1}' for i in mutations],fontsize = 165)
    plt.yticks(fontsize = 165)
    plt.grid(axis = 'x')
    plt.legend([bp['medians'][0]],['Median'],frameon=False,fontsize=150)

    plt.savefig(f'../Plots/Mutation effect on evolvability/{typ}')
    plt.close()
