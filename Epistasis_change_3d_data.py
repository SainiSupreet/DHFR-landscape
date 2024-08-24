from bio_functions import *
import scipy.io as sio

categories = ['Positive','Negative','Single Sign','Reciprocal Sign',
              'Other Sign','No Epistasis']

for bg in ['Functional','Non Functional','All']:
    for cat in categories:
        x = cat
        if cat == 'No Epistasis':
            x = 'No'
            
        data = pd.read_csv(f'../Epistasis/Epistasis Change/{x} Epistasis Data - {bg}.csv')

        ans = []
        
        for mut in range(len(categories)):
            count = get_count_from_list(data[categories[mut]].values)
            norm = sum(count.values())

            for i in range(0,22):
                if i not in count.keys():
                    count[i] = 0

            ans.append((np.array([count[i] for i in range(0,22)])/norm).tolist())

        x = ans

        sio.savemat(f'../Epistasis/Epistasis Change/3D_data_{bg}_{cat}_background.mat',dict(x=x))

for bg in ['Functional','Non Functional','All']:
    for mut in categories:
    
        ans = []

        for cat in range(len(categories)):
            s = categories[cat]
            if s == 'No Epistasis':
                s = 'No'
                
            data = pd.read_csv(f'../Epistasis/Epistasis Change/{s} Epistasis Data - {bg}.csv')
            count = get_count_from_list(data[mut])
            norm = sum(count.values())

            for i in range(0,22):
                if i not in count.keys():
                    count[i] = 0
                    
            ans.append((np.array([count[i] for i in range(0,22)])/norm).tolist())

        x = ans

        sio.savemat(f'../Epistasis/Epistasis Change/3D_data_{bg}_{mut}_mutation.mat',dict(x=x))
