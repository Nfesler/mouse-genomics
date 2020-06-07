"""Scan data set to identify the number of SNPs with allele indentification for each strain"""


import csv
import pandas as pd
from progress.bar import ShadyBar


filename_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Génome souris\\Mouse Chr {0}.txt'
out_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\Liste\\result Chr{0}.csv'



full_results = pd.DataFrame()

MIN_SAMPLES = 5  # Choix du nombre de lignée minimum données pour considérer le résultat intéressant




def analyse(filename, nb_chr):

    def line_count(filename):
        with open (filename) as file:
            file_reader = csv.reader(file)
            for i, l in enumerate(file_reader):
                pass
            return (i + 1)


    with open(filename, "r") as file:

        file_reader = csv.reader(file, delimiter='\t')
        index = 0       # numéro de la ligne
        LastSNP  = ""
        total = 0
        snps_by_strain =  dict()
        global full_results

        total_lines = line_count(filename)
        scanned_line = 0
        bar = ShadyBar(nb_chr, max= total_lines / 1000, suffix='%(percent)d%%')

        for row in file_reader:
            if index == 0:
                headers = row
                index = index + 1
            elif row[0] != LastSNP:   # Enlever les SNP écrits 2 fois de suite avec gene ID !=
                # Look for differences

                # parse = dictionnaire des lignées pour chaque allele
                n_samples = 0

                for i in range (8, len(row) - 1):
                    if row[i] != '' and not "conflict" in row[i]: # Parfois ligne avec conflit sans forcément commmencer par un ?
                        snps_by_strain[headers[i]] = snps_by_strain.get(headers[i], 0) + 1

                index = index + 1

            LastSNP = row[0]

            scanned_line = scanned_line + 1
            if scanned_line % 1000 == 0:
                bar.next()

        bar.finish()

        snps_by_strain['_MAX_'] = index
        # Use pandas library to build histogram
        pd.set_option('precision', 0)
        all_strains = pd.DataFrame(headers[8:len(row)-1])
        all_strains.columns = ['strain']

        df = pd.DataFrame.from_dict(snps_by_strain, orient="index")
        df.columns = [nb_chr]
        print(df)

#         dt = pd.DataFrame(df.groupby('strain').count()[['SNP']])
#
#         dt = pd.merge(left = all_strains, right = dt, on='strain', how='outer')
#         dt = dt.fillna(0)
#         print(dt)
#         print("Total de SNP:    ", total)
        if full_results.empty:
            full_results = df
        else:
            full_results = pd.concat([full_results, df], axis=1, sort=False)

for i in range(19):
    index = i + 1
    filename = filename_pattern.format(index)
    analyse(filename, nb_chr ='chr' + str(index))

full_results = full_results.sort_index()
print(full_results)
full_results.to_csv('C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\snps_by_strain.csv', decimal = ',')
