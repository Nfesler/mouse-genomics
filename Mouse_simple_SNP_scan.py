import csv
import pandas as pd
from progress.bar import ShadyBar



filename_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Génome souris\\Mouse Chr {0}.txt'
out_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\Liste\\result Chr{0}.csv'



full_results = pd.DataFrame()

MIN_SAMPLES = 5  # Choix du nombre de lignée minimum données pour considérer le résultat intéressant




def analyse(filename, out, nb_chr):

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
        results =  []
        total_lines = line_count(filename)
        global full_results
        scanned_line = 0
        bar = ShadyBar(nb_chr, max= total_lines / 1000, suffix='%(percent)d%%')

        for row in file_reader:
            if index == 0:
                headers = row
                index = index + 1
            elif row[0] != LastSNP:   # Enlever les SNP écrits 2 fois de suite avec gene ID !=
                # Look for differences
                parse = dict()
                # parse = dictionnaire des lignées pour chaque allele
                n_samples = 0

                for i in range (8, len(row) - 1):
                    if row[i] != '' and not "conflict" in row[i]: # Parfois ligne avec conflit sans forcément commmencer par un ?
                        n_samples = n_samples + 1
                        if row[i] in parse:
                            parse[row[i]].append(headers[i])
                        else:
                            parse[row[i]] = [ headers[i] ]

                    # Check if there is a unique element - une lignée pour un allele donné
                    onlyOnce = True
                    differenciator = ""

                    if len(parse) > 1 and n_samples > MIN_SAMPLES:
                    # On élimine les cas où on a qu'un seul type d'allele ou un nombre insuffisant de lignées pour que ce soit intéressant
                        for k,v in parse.items():
                            if len(v) == 1:
                                if differenciator == "":
                                    differenciator = k
                                else:
                                    # La deuxième fois qu'on rencontre un allele avec une seule
                                    # lignée, differenciator <> "" donc on positionne le flag onlyOnce à False
                                    # car cela veut dire qu'on a plusieurs alleles avec une seule lignée
                                    onlyOnce = False
                                    break

                if onlyOnce and (differenciator != ""):
                    result = dict()
                    result['line'] = index
                    result['SNP'] = row[0]
                    result['strain'] = parse[differenciator][0]
                    result['n_samples'] = n_samples
                    results.append(result)

                    total = total + 1

                index = index + 1

            LastSNP = row[0]

            scanned_line = scanned_line + 1
            if scanned_line % 1000 == 0:
                bar.next()

        bar.finish()

        # Use pandas library to build histogram
        pd.set_option('precision', 0)
        all_strains = pd.DataFrame(headers[8:len(row)-1])
        all_strains.columns = ['strain']
        df = pd.DataFrame(results)
        df.to_csv(out)
        dt = pd.DataFrame(df.groupby('strain').count()[['SNP']])
        dt.columns = [nb_chr]
        dt = pd.merge(left = all_strains, right = dt, on='strain', how='outer')
        dt = dt.fillna(0)
        print(dt)
        print("Total de SNP:    ", total)
        if full_results.empty:
            full_results = dt
        else:
            full_results = pd.merge(left=full_results, right=dt, on='strain', how='outer')


for i in range(19):
    index = i + 1
    filename = filename_pattern.format(index)
    out = out_pattern.format(index)
    analyse(filename, out, nb_chr ='chr' + str(index))

print(full_results)
full_results.to_csv('C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\full_results_simple.csv', decimal = ',')
