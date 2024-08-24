from libraries import *

plt.rcParams["figure.figsize"] = (50,50)
plt.rcParams['axes.linewidth'] = 20
plt.tick_params(axis='both', length = 45,width=45)
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

binn = 0.05

# Plot for Peak count - functional
path = 'Peak count - functional/'
d = []
nk = []
p = []

p_values = {
    'Landscape size':[],
    'p-value':[]
            }

for i in range(1,10):
    
    file = '../'+path+f'Peak count for {i} Point landscapes.csv'

    data1 = pd.read_csv(file)
    data2 = pd.read_csv('../'+path + f'Expected peak count for {i} Point landscapes.csv')
    
    dist_emp = []
    dist_nk = []

    for k in data1['Number of Peaks'].values:
        dist_emp += [k for j in range(data1['Number of Landscapes']\
                                      [data1['Number of Peaks']==k].values[0])]

    for k in data2['Number of Peaks'].values:
        dist_nk += [k for j in range(data2['Number of Landscapes']\
                                      [data2['Number of Peaks']==k].values[0])]


    p_values['Landscape size'].append(i)
    p_values['p-value'].append(stats.mannwhitneyu(dist_emp,dist_nk).pvalue)

    temp = pd.read_csv('../'+path+f'Detailed data/Detailed peaks in {i} \
point landscapes.csv')
    d.append(temp['Number of Peaks'].values)
    p.append(4**i)
    nk.append(temp['Predicted (NK) peaks'].values)
    
    stem1 = plt.stem([i+binn for i in data1['Number of Peaks']],
                    data1['Fraction of landscapes'],label = 'Emperical peaks',
                    linefmt='red', markerfmt='x')
    stem2 = plt.stem([i-binn for i in data2['Number of Peaks']],
                    data2['Fraction of landscapes'],label = 'NK landscape peaks',
                    linefmt='blue', markerfmt='x')
    
    plt.setp(stem1[1], 'linewidth', 15)
    stem1[0].set_markeredgewidth(90)

    plt.setp(stem2[1], 'linewidth', 15)
    stem2[0].set_markeredgewidth(90)

    plt.legend(fontsize = 150,frameon = False,loc=1)

    plt.xlabel('Number of peaks',fontsize = 175)
    plt.ylabel('Fraction of landscapes',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)

    plt.locator_params(axis='x', nbins=6)
    
    plt.savefig(f'../Plots/{path}Joint Peak count for {i} Point landscapes.png',bbox_inches='tight')
    plt.close()

df = pd.DataFrame(p_values)
df.to_csv('../Peak Probability - functional/p-values for maximally rugged landscapes.csv',index = False)
