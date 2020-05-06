import csv
import pandas as pd

limited_results = pd.DataFrame()

with open("Mouse Chr lim.txt", "r") as file:

    file_reader = csv.reader(file, delimiter='\t')
    index = 0       # numéro de la ligne
    results =  []
    LastSNP = ""
    min = 5  # Choix du nombre de lignée minimum données pour considérer le résultat intéressant
    total = 0
    for row in file_reader:
        if index == 0:
            headers = row
            index = index + 1
        elif row[0] != LastSNP:   # Enlever les SNP écrits 2 fois de suite avec gene ID !=
                                  # TO DO: Verifier que les lignes sont bien pareil??
            # Look for differences
            parse = dict()
            n_samples = 0

            for i in range (8, len(row) - 1):
                if row[i] != '' and not "conflict" in row[i]:  # row[i].startswith("?"):
                                                               # Parfois ligne avec conflit sans forcément commmencer par un ?
                    n_samples = n_samples + 1
                    if row[i] in parse:
                        parse[row[i]].append(headers[i])
                    else:
                        parse[row[i]] = [ headers[i] ]
                # Check if there is a unique element - une lignée pour un allele donné
                # parse = dictionnaire des lignées pour chaque allele
                onlyOnce = True
                differenciator = ""
                if len(parse) > 1 and n_samples > min:         # On élimine les cas où on a qu'un seul type d'allele
                                                                # ou un nombre insuffisant de lignées pour que ce soit
                                                                # intéressant
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

    pd.set_option('precision', 0)
    all_strains = pd.DataFrame(headers[8:len(row)-1])
    all_strains.columns = ['strain']
    print(all_strains)
    # Use pandas library to build histogram
    df = pd.DataFrame(results)
    df.to_csv("Result Chr lim.txt")
    dt = pd.DataFrame(df.groupby('strain').count()[['SNP']])
    dt.columns = ["Chr lim"]
    print(dt)
    print("Total de SNP:    ", total)
    dt = pd.merge(left = all_strains, right = dt, on='strain', how='outer')
    dt = dt.fillna(0)
    # dt = dt.reset_index(drop = True)
    print(dt)
    # if limited_results.empty:
    #     limited_results = dt
    # else:
    #     limited_results = pd.merge(left=limited_results, right=ft, on='lignee', how='outer')
    print("Total de SNP:    ", total)
    pd.set_option('precision', 0)
    dt.to_csv('Tableau recapchr lim.csv', decimal = ",")
