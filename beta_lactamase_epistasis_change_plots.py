from libraries import *

plt.rcParams["figure.figsize"] = (60,50)
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

data = pd.read_csv('../Beta lactamase/Epistasis change data.csv')
epistasis = ['Positive','Negative','Sign','No']

for typ in epistasis:
    
    reduced = data[data['Background Epistasis'] == typ]
    
    bp = plt.boxplot([reduced[i] for i in epistasis],
                     medianprops=dict(linestyle='-', linewidth=15, color='red'))

    plt.xlabel('Mutant Epistasis',fontsize = 175)
    plt.ylabel('Number of neighbours',fontsize = 175)
    plt.xticks(range(1,5),['Pos','Neg','Sign','No'],fontsize = 165)
    plt.yticks(range(4),range(4),fontsize = 165)
    
##    legend = plt.legend([bp['medians'][0]],
##               ['Median'],
##               frameon=True,fontsize=150)

    plt.savefig(f'../Plots/Beta lactamase/{typ} epistasis background')
    plt.close()
