from _peaks import *

path = '../Peak count - all points/'

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
    peaks,NK = get_peak_counts(i,save_data = True,
                            consider_threshold=False,
                            path=path+'Detailed data/')
    
    counts = get_count_from_list(peaks)
    total_landscapes = len(get_n_nucleotide_landscapes(i))
    
    data = pd.DataFrame({'Number of Peaks':counts.keys(),
                         'Number of Landscapes':counts.values(),
                         'Fraction of landscapes':[count/total_landscapes for count in list(counts.values())]})
    
    data = data.sort_values(by = 'Number of Peaks')
    data.to_csv(f'{path}Peak count for {i} Point landscapes.csv',
                index = False)

    NK_count = get_count_from_list(NK)

    NK_count = pd.DataFrame({'Number of Peaks':NK_count.keys(),
                         'Number of Landscapes':NK_count.values(),
                         'Fraction of landscapes':[count/total_landscapes for count in list(NK_count.values())]})
    
    NK_count = NK_count.sort_values(by = 'Number of Peaks')
    NK_count.to_csv(f'{path}Expected peak count for {i} Point landscapes.csv',
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
    'Average number of Peaks':Average_peaks,
    'Median number of Peaks':Median_peaks,
    'Variance in number of Peaks':Variance,
    'Expected average number of Peaks':NK_average,
    'Expected median number of Peaks':NK_median,
    'Expected variance number of in Peaks':NK_variance
}

compiled_df = pd.DataFrame(data = compiled_data)
compiled_df = compiled_df.to_csv(f'{path}Compiled peak data.csv',index = False)
