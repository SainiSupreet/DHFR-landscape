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

s = 7000

# Plot for Peak count - all variants
path = 'Peak count - all points/'
d = []
p = []
for i in range(1,10):
    
    file = '../'+path+f'Peak count for {i} Point landscapes.csv'
    
    data = pd.read_csv(file)

    temp = pd.read_csv('../'+path+f'Detailed data/Detailed peaks in {i} \
point landscapes.csv')
    d.append(temp['Number of Peaks'].values)
    p.append(4**i)
    
    stem = plt.stem(data['Number of Peaks'],data['Fraction of landscapes'],
                    linefmt='red', markerfmt='x')
    plt.setp(stem[1], 'linewidth', 15)
    stem[0].set_markeredgewidth(90)
    
    plt.xlabel('Number of peaks',fontsize = 175)
    plt.ylabel('Fraction of landscapes',fontsize = 175)
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)

    plt.locator_params(axis='x', nbins=6)

    plt.savefig(f'../Plots/{path}Peak count for {i} Point landscapes.png',bbox_inches='tight')
    plt.close()

bp = plt.boxplot(d,showmeans = True,meanline = True,
                 showfliers = False,
            medianprops=dict(linestyle='-', linewidth=15, color='black'),
            meanprops=dict(linestyle='--', linewidth=15, color='red'))

sct = plt.scatter(range(1,10),p,s=s,label = 'Landscape size',marker = '*',color = 'blue')

plt.xlabel('n-base pair landscape',fontsize = 175)
plt.ylabel('Number of peaks',fontsize = 175)
plt.xticks(fontsize = 165)
plt.yticks(fontsize = 165)
plt.yscale('log')

plt.legend([bp['means'][0],bp['medians'][0],sct],
           ['Mean','Median','Landscape size'],
           frameon=False,fontsize=150,loc = 'upper left')

plt.savefig('../Plots/'+path+'Compiled data')
plt.close()

# Plot for Peak count - functional
path = 'Peak count - functional/'
d = []
nk = []
p = []
for i in range(1,10):
    
    file = '../'+path+f'Peak count for {i} Point landscapes.csv'

    data = pd.read_csv(file)

    temp = pd.read_csv('../'+path+f'Detailed data/Detailed peaks in {i} \
point landscapes.csv')
    d.append(temp['Number of Peaks'].values)
    p.append(4**i)
    nk.append(temp['Predicted (NK) peaks'].values)
    
    stem = plt.stem(data['Number of Peaks'],data['Fraction of landscapes'],
                    linefmt='red', markerfmt='x')
    plt.setp(stem[1], 'linewidth', 15)
    stem[0].set_markeredgewidth(90)

    plt.xlabel('Number of peaks',fontsize = 175)
    plt.ylabel('Fraction of landscapes',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)

    plt.locator_params(axis='x', nbins=6)
    
    plt.savefig(f'../Plots/{path}Peak count for {i} Point landscapes.png',bbox_inches='tight')
    plt.close()

    data = pd.read_csv('../'+path + f'Expected peak count for {i} Point landscapes.csv')

    stem = plt.stem(data['Number of Peaks'],data['Fraction of landscapes'],
                    linefmt='blue', markerfmt='x')
    plt.setp(stem[1], 'linewidth', 15)
    stem[0].set_markeredgewidth(90)
    
    plt.xlabel('Number of peaks',fontsize = 175)
    plt.ylabel('Fraction of landscapes',fontsize = 175)

    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)

    plt.locator_params(axis='x', nbins=6)

    plt.savefig(f'../Plots/{path}Predicted Peak count for {i} Point landscapes.png',bbox_inches='tight')
    plt.close()

bp1 = plt.boxplot(d,showmeans = True,meanline = True,showfliers = False,
            medianprops=dict(linestyle='-', linewidth=15, color='black'),
            meanprops=dict(linestyle='--', linewidth=15, color='red'))

bp2 = plt.boxplot(nk,
            medianprops=dict(linestyle=':', linewidth=15, color='blue'),
            showbox=False, showcaps=False,showfliers = False,whis = 0)

sct = plt.scatter(range(1,10),p,s=s,label = 'Landscape size',marker = '*',color = 'blue')

plt.xlabel('n-base pair landscape',fontsize = 175)
plt.ylabel('Number of peaks',fontsize = 175)
plt.xticks(fontsize = 165)
plt.yticks(fontsize = 165)
plt.yscale('log')

plt.legend([bp1['means'][0],bp1['medians'][0],bp2['medians'][0],sct],
           ['Mean','Median','NK model median','Landscape size'],
           frameon=False,fontsize=150,loc = 2)

plt.savefig('../Plots/'+path+'Compiled data')
plt.close()
    
# Plot for Peak Probability - all points
path = 'Peak Probability - all points/'
d = []
nk = []
for i in range(1,10):
    
    file = '../'+path+f'Peak Probability for {i} Point landscapes.csv'

    data = pd.read_csv(file)

    temp = pd.read_csv('../'+path+f'Detailed data/Detailed peak probability in {i} \
point landscapes.csv')
    d.append(temp['Peak Probability'].values)

    stem = plt.stem(data['Peak Probability'],data['Fraction of landscapes'],
                    linefmt='red', markerfmt='x')
    plt.setp(stem[1], 'linewidth', 15)
    stem[0].set_markeredgewidth(90)
    
    plt.xlabel('Peak probability',fontsize = 175)
    plt.ylabel('Fraction of landscapes',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)

    plt.locator_params(axis='x', nbins=6)
    
    plt.savefig(f'../Plots/{path}Peak Probability for {i} Point landscapes.png',bbox_inches='tight')
    plt.close()

bp = plt.boxplot(d,showmeans = True,meanline = True,showfliers = False,
            medianprops=dict(linestyle='-', linewidth=15, color='black'),
            meanprops=dict(linestyle='--', linewidth=15, color='red'))

plt.legend([bp['means'][0],bp['medians'][0]],
           ['Mean','Median'],
           frameon=False,fontsize=150,
           loc = 1)
plt.xlabel('n-base pair landscape',fontsize = 175)
plt.ylabel('Peak probability',fontsize = 175)

plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.yscale('log')

plt.savefig('../Plots/'+path+'Compiled data')
plt.close()

# Plot for Peak Probability - functional
path = 'Peak Probability - functional/'
d = []
nk = []
for i in range(1,10):
    
    file = '../'+path+f'Peak Probability for {i} Point landscapes.csv'

    data = pd.read_csv(file)

    temp = pd.read_csv('../'+path+f'Detailed data/Detailed peak probability in {i} \
point landscapes.csv')
    d.append(temp['Peak Probability'].values)
    nk.append(temp['Expected (NK) probability'].values)

    stem = plt.stem(data['Peak Probability'],data['Fraction of landscapes'],
                    linefmt='red', markerfmt='x')
    plt.setp(stem[1], 'linewidth', 15)
    stem[0].set_markeredgewidth(90)
    
    plt.xlabel('Peak probability',fontsize = 175)
    plt.ylabel('Fraction of landscapes',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)

    plt.locator_params(axis='x', nbins=6)
    
    plt.savefig(f'../Plots/{path}Peak Probability for {i} Point landscapes.png',bbox_inches='tight')
    plt.close()

    data = pd.read_csv('../'+path + f'Expected peak probability for {i} Point landscapes.csv')

    stem = plt.stem(data['Peak Probability'],data['Fraction of landscapes'],
                    linefmt='blue', markerfmt='x')
    plt.setp(stem[1], 'linewidth', 15)
    stem[0].set_markeredgewidth(90)
    
    plt.xlabel('Expected peak probability',fontsize = 175)
    plt.ylabel('Fraction of landscapes',fontsize = 175)

    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)

    plt.locator_params(axis='x', nbins=6)
    
    plt.savefig(f'../Plots/{path}Predicted Peak Probability for {i} Point landscapes.png',bbox_inches='tight')
    plt.close()

bp1 = plt.boxplot(d,showmeans = True,meanline = True,showfliers = False,
            medianprops=dict(linestyle='-', linewidth=15, color='black'),
            meanprops=dict(linestyle='--', linewidth=15, color='red'))

bp2 = plt.boxplot(nk,
            medianprops=dict(linestyle=':', linewidth=15, color='blue'),
            showbox=False, showcaps=False,showfliers = False,whis = 0)

plt.legend([bp1['means'][0],bp1['medians'][0],bp2['medians'][0]],
           ['Mean','Median','NK model median'],
           frameon=False,fontsize=150,loc = 1)
plt.xlabel('n-base pair landscape',fontsize = 175)
plt.ylabel('Peak probability',fontsize = 175)

plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.yscale('log')

plt.savefig('../Plots/'+path+'Compiled data')
plt.close()
