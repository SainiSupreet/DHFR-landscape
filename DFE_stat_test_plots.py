from _DFE import *
from _mutation_effect import *

plt.rcParams["figure.figsize"] = (100,50)
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

path = '../DFE/'
background_bins = dfe_bins(list(fitness_data.values()),list(fitness_data.keys()),bin_size=0.05)
bin_values = list(background_bins.keys())
bin_values.sort()

facecolor = '#E8E8E8'
shift = -0

# KS Test
data = []
for binV in bin_values:
    try:
        df = pd.read_csv(f'../DFE/KS test/Background Fitness {binV}.csv')
        data.append(df['p value'].values)
    except FileNotFoundError:
        continue

##bp = plt.boxplot(data,showmeans = True,
##                 meanline = True,showfliers = False,patch_artist=True,
##                 medianprops=dict(linestyle='-', linewidth=15, color='black'),
##                 meanprops=dict(linestyle='-', linewidth=15, color='red'))
##
##[i.set_facecolor(facecolor) for i in bp['boxes']]

line = plt.plot([1,51],[0.05,0.05],
                linestyle = (0,(5,5)),color = 'blue')

violin = plt.violinplot(data,positions = np.array(list(range(1,52))) + shift,
                        side = 'low',
                        showmedians = True,
                        showextrema = False,
                        showmeans = True,
                        widths = 1)

for pc in violin['bodies']:
    #pc.set_facecolor('#D43F3A')
    pc.set_facecolor(facecolor)
    pc.set_edgecolor('black')
    pc.set_alpha(1)
    pc.set_linewidths(5)
    #pc.set_hatch('x')

violin['cmedians'].set(color = 'black')
violin['cmeans'].set(color = 'red')

##x = 1
##for dt in data:
##    plt.scatter([x for i in data[x-1]],data[x-1],color = 'darkblue')
##    x+=1

plt.xticks(range(1,52),[round(i,2) for i in bin_values[:-1]])

plt.legend([violin['cmeans'],violin['cmedians'],line[0]],
           ['Mean','Median','Critical p-value'],
           frameon=False,fontsize=150,loc = 'upper center',
           bbox_to_anchor=(0.5,1.12),ncol = 3)

plt.xlabel('Background Fitness',fontsize = 175)
plt.ylabel('KS test p value',fontsize = 175)
plt.ylim(0,1)
plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.locator_params('x',nbins = 8)
plt.tick_params(axis='both', length = 45,width=25)
plt.margins(x=0.01)

plt.savefig('../Plots/DFE/KS test p values')
plt.close()

# MWU test
data = []
for binV in bin_values:
    try:
        df = pd.read_csv(f'../DFE/MWU test/Background Fitness {binV}.csv')
        data.append(df['p value'].values)
    except FileNotFoundError:
        continue

##bp = plt.boxplot(data,showmeans = True,
##                 meanline = True,showfliers = False,patch_artist=True,
##                 medianprops=dict(linestyle='-', linewidth=15, color='black'),
##                 meanprops=dict(linestyle='-', linewidth=15, color='red'))
##
##[i.set_facecolor(facecolor) for i in bp['boxes']]

line = plt.plot([1,51],[0.05,0.05],
                linestyle = (0,(5,5)),color = 'blue')

violin = plt.violinplot(data,positions = np.array(list(range(1,52))) + shift,
                        side = 'low',
                        showmedians = True,
                        showextrema = False,
                        showmeans = True,
                        widths = 1)

for pc in violin['bodies']:
    #pc.set_facecolor('#D43F3A')
    pc.set_facecolor(facecolor)
    pc.set_edgecolor('black')
    pc.set_alpha(1)
    pc.set_linewidths(5)
    #pc.set_hatch('x')

violin['cmedians'].set(color = 'black')
violin['cmeans'].set(color = 'red')

##x = 1
##for dt in data:
##    plt.scatter([x for i in data[x-1]],data[x-1],color = 'darkblue')
##    x+=1

plt.xticks(range(1,52),[round(i,2) for i in bin_values[:-1]])

plt.legend([violin['cmeans'],violin['cmedians'],line[0]],
           ['Mean','Median','Critical p-value'],
           frameon=False,fontsize=150,loc = 'upper center',
           bbox_to_anchor=(0.5,1.12),ncol = 3)

plt.xlabel('Background Fitness',fontsize = 175)
plt.ylabel('MWU test p value',fontsize = 175)
plt.ylim(0,1)
plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.locator_params('x',nbins = 8)
plt.tick_params(axis='both', length = 45,width=25)
plt.margins(x=0.01)

plt.savefig('../Plots/DFE/MWU test p values')
plt.close()
