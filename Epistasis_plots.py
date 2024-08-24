from bio_functions import *

plt.rcParams["figure.figsize"] = (50,50)
plt.rcParams['axes.linewidth'] = 20
plt.rcParams['lines.linewidth'] = 15

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

path = '../Epistasis/'

# Fraction of Epistasis - functional
file = path+'Epistasis Fractions - Functional Variants.csv'
data = pd.read_csv(file)

bp = plt.boxplot([data['Positive Epistasis Fraction'],
                  data['Negative Epistasis Fraction'],
                  data['Reciprocal Sign Epistasis Fraction'],
                  data['Single Sign Epistasis Fraction'],
                  data['Other Sign Epistasis Fraction'],
                  data['No Epistasis Fraction']],
                 showmeans = True,meanline = True,showfliers = False,
                 medianprops=dict(linestyle='-', linewidth=15, color='black'),
                 meanprops=dict(linestyle='--', linewidth=15, color='red'))

plt.legend([bp['means'][0],bp['medians'][0]],
           ['Mean','Median'],
           frameon=False,fontsize=150)

plt.xticks(range(1,7),['PE', 'NE','RSE',
                       'SSE','OSE', 'No'],
           fontsize = 150)
plt.yticks(fontsize = 165)
plt.tick_params(axis='both', length = 45,width=10)
plt.ylim(0,1)

plt.xlabel('Nature of epistasis',fontsize = 175,labelpad = 60)
plt.ylabel('Fraction of genotypes',fontsize = 175,labelpad = 60)

plt.grid()

plt.savefig('../Plots/Epistasis/Epistasis in Functional Variants')
plt.close()

# Fraction of Epistasis - non functional
file = path+'Epistasis Fractions - Non Functional Variants.csv'
data = pd.read_csv(file)

bp = plt.boxplot([data['Positive Epistasis Fraction'],
                  data['Negative Epistasis Fraction'],
                  data['Reciprocal Sign Epistasis Fraction'],
                  data['Single Sign Epistasis Fraction'],
                  data['Other Sign Epistasis Fraction'],
                  data['No Epistasis Fraction']],
                 showmeans = True,meanline = True,showfliers = False,
                 medianprops=dict(linestyle='-', linewidth=15, color='black'),
                 meanprops=dict(linestyle='--', linewidth=15, color='red'))

plt.legend([bp['means'][0],bp['medians'][0]],
           ['Mean','Median'],
           frameon=False,fontsize=150)

plt.xticks(range(1,7),['PE', 'NE','RSE',
                       'SSE','OSE', 'No'],
           fontsize = 150)
plt.yticks(fontsize = 165)
plt.tick_params(axis='both', length = 45,width=10)
plt.ylim(0,1)

plt.xlabel('Nature of epistasis',fontsize = 175,labelpad = 60)
plt.ylabel('Fraction of genotypes',fontsize = 175,labelpad = 60)

plt.grid()

plt.savefig('../Plots/Epistasis/Epistasis in Non Functional Variants')
plt.close()

# Fraction of Epistasis - all
file = path+'Epistasis Fractions - All Variants.csv'
data = pd.read_csv(file)

bp = plt.boxplot([data['Positive Epistasis Fraction'],
                  data['Negative Epistasis Fraction'],
                  data['Reciprocal Sign Epistasis Fraction'],
                  data['Single Sign Epistasis Fraction'],
                  data['Other Sign Epistasis Fraction'],
                  data['No Epistasis Fraction']],
                 showmeans = True,meanline = True,showfliers = False,
                 medianprops=dict(linestyle='-', linewidth=15, color='black'),
                 meanprops=dict(linestyle='--', linewidth=15, color='red'))

plt.legend([bp['means'][0],bp['medians'][0]],
           ['Mean','Median'],
           frameon=False,fontsize=150)

plt.xticks(range(1,7),['PE', 'NE','RSE',
                       'SSE','OSE', 'No'],
           fontsize = 150)
plt.yticks(fontsize = 165)
plt.tick_params(axis='both', length = 45,width=10)
plt.ylim(0,1)

plt.xlabel('Nature of epistasis',fontsize = 175,labelpad = 60)
plt.ylabel('Fraction of genotypes',fontsize = 175,labelpad = 60)

plt.grid()

plt.savefig('../Plots/Epistasis/Epistasis in All Variants')
plt.close()
