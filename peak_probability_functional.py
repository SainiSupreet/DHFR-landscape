from _peaks import *

path = '../Peak Probability - functional/'

Points = []
Average_peaks = []
Median_peaks = []
Variance = []

NK_average = []
NK_median = []
NK_variance = []

print('Starting')

for i in range(1,10):
    print(f'Analysing {i} point landscapes')
    peaks,NK = get_peak_probability(i,save_data = True,
                            consider_threshold=True,
                            path=path+'Detailed data/')
    
    counts = get_count_from_list(peaks)
    total_landscapes = len(get_n_nucleotide_landscapes(i))
    
    data = pd.DataFrame({'Peak Probability':counts.keys(),
                         'Number of Landscapes':counts.values(),
                         'Fraction of landscapes':[count/total_landscapes for count in list(counts.values())]})
    
    data = data.sort_values(by = 'Peak Probability')
    data.to_csv(f'{path}Peak probability for {i} Point landscapes.csv',
                index = False)

    NK_count = get_count_from_list(NK)

    NK_count = pd.DataFrame({'Peak Probability':NK_count.keys(),
                         'Number of Landscapes':NK_count.values(),
                         'Fraction of landscapes':[count/total_landscapes for count in list(NK_count.values())]})
    
    NK_count = NK_count.sort_values(by = 'Peak Probability')
    NK_count.to_csv(f'{path}Expected peak probability for {i} Point landscapes.csv',
                    index = False)

    Points.append(i)
    Average_peaks.append(mean(peaks))
    Median_peaks.append(median(peaks))
    Variance.append(variance(peaks))

    NK_average.append(mean(NK))
    NK_median.append(median(NK))
    NK_variance.append(variance(NK))
    print('\n')

print('Saving Compiled data')

compiled_data = {
    'n-points':Points,
    'Average Peak Probability':Average_peaks,
    'Median Peak Probability':Median_peaks,
    'Variance in Peak Probability':Variance,
    'Expected average Peak Probability':NK_average,
    'Expected median Peak Probability':NK_median,
    'Expected variance Peak Probability':NK_variance
}

compiled_df = pd.DataFrame(data = compiled_data)
compiled_df = compiled_df.to_csv(f'{path}Compiled peak data.csv',index = False)
