from bio_functions import *
from _DFE import *

plt.rcParams["figure.figsize"] = (60,50)
plt.rcParams['axes.linewidth'] = 20
plt.rcParams['lines.linewidth'] = 15

plt.rcParams['xtick.major.size']=50
plt.rcParams['xtick.major.width']=10
plt.rcParams['ytick.major.size']=50
plt.rcParams['ytick.major.width']=10

plt.rcParams['grid.linewidth']=3

plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.constrained_layout.h_pad']=0.5
plt.rcParams['figure.constrained_layout.w_pad']=0.5
plt.rcParams['savefig.pad_inches'] = 0.5

colors = ['magenta','blue','red','green']
fitness = [-0.47,0.03,0.53,1.03]
alpha = [0.5,0.5,0.5]
breaks = 200

col = 0
for binn in fitness:
    
    col += 1
    
    data = pd.read_csv(f'../DFE/Probability Distributions - all/BF {binn}.csv')

    cubic_model = interpolate.interp1d(data['SC'],data['PDF'],kind = "cubic")

    values = np.linspace(min(data['SC']),max(data['SC']),breaks)

    plt.plot(values,cubic_model(values),linewidth = 15,
             label = f'BF = {round(binn,2)}',
             color = colors[(col-1)%len(colors)])

plt.plot([0,0],[0,3],linestyle = '--',color = 'black',alpha = 0.5,label = 'Neutral\nmutations')

plt.legend(frameon = False,fontsize = 150)
plt.xticks(fontsize = 165)
plt.locator_params('x',nbins = 6)
plt.yticks(fontsize = 150,rotation = 0)
plt.xlabel('Selection Coefficient',fontsize = 175)
plt.ylabel('Probability density',fontsize = 175,labelpad=60)

plt.savefig('../Plots/DFE/DFE probability distribution 2D')
