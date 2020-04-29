import pandas as pd

filename_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\StrainTab\\result 2 {0}.txt'

full_results = pd.DataFrame()

def analyse(filename, columnName):
    global full_results
    data = pd.read_table(filename, sep='\s+', skiprows=2, skipfooter=1, engine='python', header=None)
    # data = data.drop(1, axis=1)
    data.columns=['lignee', columnName]
    #print(data.columns)
    print(data)
    if full_results.empty:
        full_results = data
    else:
        full_results = pd.merge(left=full_results, right=data, on='lignee', how='outer')
    #print(full_results)

for i in range (19):
    index = i + 1
    filename = filename_pattern.format(index)
    analyse(filename, columnName='chr'+str(index))

print(full_results)
full_results.to_csv('C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\full_results.csv')
