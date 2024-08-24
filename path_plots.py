from libraries import *

target = 'AAGGAAATG'

def factorial(n):
    '''
    Returns n!
    '''
    if n > 1:
        return n*factorial(n-1)
    elif n == 0 or n == 1:
        return 1
    else:
        return 0

plt.rcParams["figure.figsize"] = (50,50)
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

categories = ['Functional','Non Functional','All']

arrowprops = dict(arrowstyle="-|>",
                  linewidth = 10,
                  color = 'black')

bbox = dict(boxstyle="Square",
            color = 'whitesmoke')

for category in categories:
    for i in range(1,10):
        plt.close()
        try:
            data = pd.read_csv(f'../Peak Accessibility/Count data for {i} distance from {target} - {category}.csv')
        except FileNotFoundError:
            break
        
        stem = plt.stem(data['Number of Accessible Paths'],
                        data['PDF'],
                        linefmt='red', markerfmt='x')

        stem[1].set_linewidth(15)
        stem[0].set_markeredgewidth(90)

        plt.xlabel('Number of accessible paths',fontsize = 175)
        plt.ylabel('Fraction of variants',fontsize = 175)
        plt.xticks(fontsize = 165)
        plt.yticks(fontsize = 165)

        plt.locator_params(axis='both', nbins=6)

        if i > 3:
            if i == 4 and category == 'Non Functional':
                assert True
            else:
                fractions = data['PDF'].values
                try:
                    max_val = max(fractions[1:])
                except:
                    break
                plt.ylim(0,max_val*1.1)
                plt.annotate(f'{round(fractions[0],2)}',xy=(0,max_val*1.09),
                             fontsize = 145,bbox=bbox, arrowprops=arrowprops,
                             xytext=(max(data['Number of Accessible Paths'].values)/6,
                                     max_val*1.038))

        plt.savefig(f'../Plots/Peak Accessibility/{category} variants - {i} distance from {target}')
        plt.close()
