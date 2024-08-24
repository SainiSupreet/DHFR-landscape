from _mutation_effect import *

plt.rcParams["figure.figsize"] = (60,50)
plt.rcParams['axes.linewidth'] = 20
plt.rcParams['lines.linewidth'] = 10
plt.rcParams['xtick.major.size']=50
plt.rcParams['xtick.major.width']=10
plt.rcParams['ytick.major.size']=50
plt.rcParams['ytick.major.width']=10

plt.rcParams['axes.labelpad'] = 40

plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.constrained_layout.h_pad']=0.5
plt.rcParams['figure.constrained_layout.w_pad']=0.5
plt.rcParams['savefig.pad_inches'] = 0.5

# Plot all
path = 'Mutation effect - all points/'

all_points = {'BF':[],'SC':[]}
nf = {'BF':[],'SC':[]}
fun = {'BF':[],'SC':[]}

lin_reg = pd.read_csv(f'../{path}Linear regression.csv')

for mutations in generate_mutations(1):
    mutation = mutations[0]

    all_data = pd.read_csv(f'../{path}all/{mutation[0]} at {mutation[1]+1} mutation data.csv')
    fun_data = pd.read_csv(f'../{path}functional/{mutation[0]} at {mutation[1]+1} mutation data.csv')
    nf_data = pd.read_csv(f'../{path}non functional/{mutation[0]} at {mutation[1]+1} mutation data.csv')

    all_points['BF']+=list(all_data['Background Fitness'].values)
    all_points['SC']+=list(all_data['SC'].values)

    fun['BF']+=list(fun_data['Background Fitness'].values)
    fun['SC']+=list(fun_data['SC'].values)

    nf['BF']+=list(nf_data['Background Fitness'].values)
    nf['SC']+=list(nf_data['SC'].values)


    # Functional
    plt.hexbin(fun_data['Background Fitness'],fun_data['SC'],
               bins = 'log',gridsize = 30)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
    cbar.set_label('Density',fontsize = 175)
    
    plt.plot([max(fun_data['Background Fitness']),
              min(fun_data['Background Fitness'])],
             [0,0],'--',color = 'black')
    
    plt.locator_params(axis='x', nbins=6)
    plt.xlabel('Background Fitness',fontsize = 175)
    plt.ylabel('Selection Coefficient',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)
    plt.savefig(f'../Plots/{path}Functional {mutation[0]} at {mutation[1]+1}')
    plt.close()

    # Non functional
    plt.hexbin(nf_data['Background Fitness'],nf_data['SC'],
               bins = 'log',gridsize = 30)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
    cbar.set_label('Density',fontsize = 175)
    
    plt.plot([max(nf_data['Background Fitness']),
              min(nf_data['Background Fitness'])],
             [0,0],'--',color = 'black')
    
    plt.locator_params(axis='x', nbins=6)
    plt.xlabel('Background Fitness',fontsize = 175)
    plt.ylabel('Selection Coefficient',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)
    plt.savefig(f'../Plots/{path}Non functional {mutation[0]} at {mutation[1]+1}')
    plt.close()

    # All
    slope = lin_reg[lin_reg['Mutation base'] == mutation[0]]\
            [lin_reg['Mutation locus'] == mutation[1]+1]\
            ['Linear slope'].values
    
    intercept = lin_reg[lin_reg['Mutation base'] == mutation[0]]\
                [lin_reg['Mutation locus'] == mutation[1]+1]\
                ['Intercept'].values
    
    plt.hexbin(all_data['Background Fitness'],all_data['SC'],
               bins = 'log',gridsize = 30)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
    cbar.set_label('Density',fontsize = 175)
    
    plt.plot([max(all_data['Background Fitness']),
              min(all_data['Background Fitness'])],
             [0,0],'--',color = 'black',linewidth=15)

    plt.scatter([-intercept[0]/slope[0]],[0],s=1000,color = 'red')

    x_range = np.array([max(all_data['Background Fitness']),min(all_data['Background Fitness'])])
    plt.plot(x_range,
             x_range*slope+intercept,'-',color = 'red')

    plt.locator_params(axis='x', nbins=6)
    plt.xlabel('Background Fitness',fontsize = 175)
    plt.ylabel('Selection Coefficient',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)
    plt.savefig(f'../Plots/{path}All {mutation[0]} at {mutation[1]+1}')
    plt.close()

# Functional

slope = lin_reg[lin_reg['Mutation base'] == 'Functional']\
        [lin_reg['Mutation locus'] == 'All']\
        ['Linear slope'].values

intercept = lin_reg[lin_reg['Mutation base'] == 'Functional']\
            [lin_reg['Mutation locus'] == 'All']\
            ['Intercept'].values

plt.hexbin(fun['BF'],fun['SC'],
           bins = 'log',gridsize = 30)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
cbar.set_label('Density',fontsize = 175)

plt.plot([max(fun['BF']),
          min(fun['BF'])],
         [0,0],'--',color = 'black')

plt.scatter([-intercept[0]/slope[0]],[0],s=1000,color = 'red')

x_range = np.array([max(fun['BF']),min(fun['BF'])])
plt.plot(x_range,
         x_range*slope+intercept,'-',color = 'red')

plt.locator_params(axis='x', nbins=6)
plt.xlabel('Background Fitness',fontsize = 175)
plt.ylabel('Selection Coefficient',fontsize = 175)

plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.savefig(f'../Plots/{path}Compiled Functional')
plt.close()

# Non functional

slope = lin_reg[lin_reg['Mutation base'] == 'Non functional']\
        [lin_reg['Mutation locus'] == 'All']\
        ['Linear slope'].values

intercept = lin_reg[lin_reg['Mutation base'] == 'Non functional']\
            [lin_reg['Mutation locus'] == 'All']\
            ['Intercept'].values

plt.hexbin(nf['BF'],nf['SC'],
           bins = 'log',gridsize = 30)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
cbar.set_label('Density',fontsize = 175)

plt.plot([max(nf['BF']),
          min(nf['BF'])],
         [0,0],'--',color = 'black')

plt.scatter([-intercept[0]/slope[0]],[0],s=1000,color = 'red')

x_range = np.array([max(nf['BF']),min(nf['BF'])])
plt.plot(x_range,
         x_range*slope+intercept,'-',color = 'red')

plt.locator_params(axis='x', nbins=6)
plt.xlabel('Background Fitness',fontsize = 175)
plt.ylabel('Selection Coefficient',fontsize = 175)

plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.savefig(f'../Plots/{path}Compiled non functional')
plt.close()


# All
plt.hexbin(all_points['BF'],all_points['SC'],
           bins = 'log',gridsize = 30)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
cbar.set_label('Density',fontsize = 175)

plt.plot([max(all_points['BF']),
          min(all_points['BF'])],
         [0,0],'--',color = 'black')

plt.locator_params(axis='x', nbins=6)
plt.xlabel('Background Fitness',fontsize = 175)
plt.ylabel('Selection Coefficient',fontsize = 175)

plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.savefig(f'../Plots/{path}Compiled all')
plt.close()

# Pivot point
plt.figure(figsize = (325,70))
data = lin_reg['Pivot point (background fitness)'][:-2]
data = data.values

#Removing the outlier
normal_vals = data[:6].tolist()+data[7:].tolist()
reduced_data = np.array(data[:6].tolist() + [mean(normal_vals)] + data[7:].tolist())

plt.bar([f'{i[0]} -> {i[1]} at {i[2]+1}' for i in expanded_mutations],
         reduced_data-mean(reduced_data),label = 'Pivot points',bottom = mean(reduced_data),color = 'teal')

outlier = plt.bar('A -> G at 9',data[6] - mean(reduced_data),label = 'Outlier',bottom = mean(reduced_data),
                   hatch = '/',fill = False,edgecolor = 'black',linewidth = 10,linestyle = '--')

plt.ylim(-2.5,0.5)

plt.plot([-1,len(data)],[threshold,threshold],'--',color = 'red',
         linewidth = 20,label = 'Cutoff fitness')
plt.plot([-1,len(data)],[mean(reduced_data),mean(reduced_data)],'-',
         color = 'black',linewidth = 20,label='Mean pivot point')

plt.xticks(fontsize=165,rotation = 90)
plt.yticks(fontsize=165)
plt.ylabel('Background Fitness',fontsize = 175)
plt.xlabel('Mutation',fontsize = 175)
plt.locator_params(axis='y', nbins=6)
plt.legend(fontsize = 150,frameon = False,loc=0)
plt.margins(x=0.005)
plt.savefig(f'../Plots/Linear Regression')

for mutation in expanded_mutations:
    ind_data = pd.read_csv(f'../Mutation effect - all points/Expanded mutations - all (for linear regression)/{mutation[0]}->{mutation[1]} at {mutation[2]+1}.csv')
    
    slope = lin_reg[lin_reg['Mutation base'] == f'{mutation[0]} -> {mutation[1]}']\
            [lin_reg['Mutation locus'] == mutation[2]+1]\
            ['Linear slope'].values
    
    intercept = lin_reg[lin_reg['Mutation base'] == f'{mutation[0]} -> {mutation[1]}']\
            [lin_reg['Mutation locus'] == mutation[2]+1]\
                ['Intercept'].values
    
    plt.hexbin(ind_data['Background Fitness'],ind_data['SC'],
               bins = 'log',gridsize = 30)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
    cbar.set_label('Density',fontsize = 175)
    
    plt.plot([max(ind_data['Background Fitness']),
              min(ind_data['Background Fitness'])],
             [0,0],'--',color = 'black',linewidth=15)

    plt.scatter([-intercept[0]/slope[0]],[0],s=1000,color = 'red')

    x_range = np.array([min(ind_data['Background Fitness']),
                        max(ind_data['Background Fitness'])])
    plt.plot(x_range,
             x_range*slope+intercept,'-',color = 'red')

    plt.locator_params(axis='x', nbins=6)
    plt.xlabel('Background Fitness',fontsize = 175)
    plt.ylabel('Selection Coefficient',fontsize = 175)
    plt.xlim(x_range)
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)
    plt.savefig(f'../Plots/Mutation effect - all points/Linear regression/{mutation[0]} -> {mutation[1]} at {mutation[2]+1}')
    plt.close()

# Plot only functional
path = 'Mutation effect - functional/'

all_points = {'BF':[],'SC':[]}
nf = {'BF':[],'SC':[]}
fun = {'BF':[],'SC':[]}

lin_reg = pd.read_csv(f'../{path}Linear regression.csv')

for mutations in generate_mutations(1):
    mutation = mutations[0]

    all_data = pd.read_csv(f'../{path}all/{mutation[0]} at {mutation[1]+1} mutation data.csv')
    fun_data = pd.read_csv(f'../{path}functional/{mutation[0]} at {mutation[1]+1} mutation data.csv')
    nf_data = pd.read_csv(f'../{path}non functional/{mutation[0]} at {mutation[1]+1} mutation data.csv')

    all_points['BF']+=list(all_data['Background Fitness'].values)
    all_points['SC']+=list(all_data['SC'].values)

    fun['BF']+=list(fun_data['Background Fitness'].values)
    fun['SC']+=list(fun_data['SC'].values)

    nf['BF']+=list(nf_data['Background Fitness'].values)
    nf['SC']+=list(nf_data['SC'].values)

    # Functional
    plt.hexbin(fun_data['Background Fitness'],fun_data['SC'],
               bins = 'log',gridsize = 30)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
    cbar.set_label('Density',fontsize = 175)
    
    plt.plot([max(fun_data['Background Fitness']),
              min(fun_data['Background Fitness'])],
             [0,0],'--',color = 'black')

    plt.locator_params(axis='x', nbins=6)
    plt.xlabel('Background Fitness',fontsize = 175)
    plt.ylabel('Selection Coefficient',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)
    plt.savefig(f'../Plots/{path}Functional {mutation[0]} at {mutation[1]+1}')
    plt.close()

    # Non functional
    plt.hexbin(nf_data['Background Fitness'],nf_data['SC'],
               bins = 'log',gridsize = 30)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
    cbar.set_label('Density',fontsize = 175)
    
    plt.plot([max(nf_data['Background Fitness']),
              min(nf_data['Background Fitness'])],
             [0,0],'--',color = 'black')

    plt.locator_params(axis='x', nbins=6)
    plt.xlabel('Background Fitness',fontsize = 175)
    plt.ylabel('Selection Coefficient',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)
    plt.savefig(f'../Plots/{path}Non functional {mutation[0]} at {mutation[1]+1}')
    plt.close()

    # All
    slope = lin_reg[lin_reg['Mutation base'] == mutation[0]][lin_reg['Mutation locus'] == mutation[1]+1]\
            ['Linear slope'].values
    intercept = lin_reg[lin_reg['Mutation base'] == mutation[0]][lin_reg['Mutation locus'] == mutation[1]+1]\
            ['Intercept'].values
    
    plt.hexbin(all_data['Background Fitness'],all_data['SC'],
               bins = 'log',gridsize = 30)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
    cbar.set_label('Density',fontsize = 175)
    
    plt.plot([max(all_data['Background Fitness']),
              min(all_data['Background Fitness'])],
             [0,0],'--',color = 'black',linewidth=15)

    plt.scatter([-intercept[0]/slope[0]],[0],s=1000,color = 'red')

    x_range = np.array([max(all_data['Background Fitness']),min(all_data['Background Fitness'])])
    plt.plot(x_range,
             x_range*slope+intercept,'-',color = 'red')

    plt.locator_params(axis='x', nbins=6)
    plt.xlabel('Background Fitness',fontsize = 175)
    plt.ylabel('Selection Coefficient',fontsize = 175)
    
    plt.xticks(fontsize=165)
    plt.yticks(fontsize=165)
    plt.savefig(f'../Plots/{path}All {mutation[0]} at {mutation[1]+1}')
    plt.close()

# Functional
plt.hexbin(fun['BF'],fun['SC'],
           bins = 'log',gridsize = 30)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
cbar.set_label('Density',fontsize = 175)

plt.plot([max(fun['BF']),
          min(fun['BF'])],
         [0,0],'--',color = 'black')

plt.locator_params(axis='x', nbins=6)
plt.xlabel('Background Fitness',fontsize = 175)
plt.ylabel('Selection Coefficient',fontsize = 175)

plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.savefig(f'../Plots/{path}Compiled Functional')
plt.close()

# Non functional
plt.hexbin(nf['BF'],nf['SC'],
           bins = 'log',gridsize = 30)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
cbar.set_label('Density',fontsize = 175)

plt.plot([max(nf['BF']),
          min(nf['BF'])],
         [0,0],'--',color = 'black')

plt.locator_params(axis='x', nbins=6)
plt.xlabel('Background Fitness',fontsize = 175)
plt.ylabel('Selection Coefficient',fontsize = 175)

plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.savefig(f'../Plots/{path}Compiled non functional')
plt.close()

# All
plt.hexbin(all_points['BF'],all_points['SC'],
           bins = 'log',gridsize = 30)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = 150,size = 30,width = 10)
cbar.set_label('Density',fontsize = 175)

plt.plot([max(all_points['BF']),
          min(all_points['BF'])],
         [0,0],'--',color = 'black')

plt.locator_params(axis='x', nbins=6)
plt.xlabel('Background Fitness',fontsize = 175)
plt.ylabel('Selection Coefficient',fontsize = 175)

plt.xticks(fontsize=165)
plt.yticks(fontsize=165)
plt.savefig(f'../Plots/{path}Compiled all')
plt.close()
