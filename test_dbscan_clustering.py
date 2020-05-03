import csv
import pandas as pd

limited_results = pd.DataFrame()

with open("Mouse Chr lim.txt", "r") as file:
    rst = open("Result Chr lim.txt" , "w")
    rst.write("\t".join(["ligne", "SNP", "allele", "lignee", "Nbre de lignee"]))
    rst.write("\n")
    file_reader = csv.reader(file, delimiter='\t')
    index = 0       # numéro de la ligne
    results =  []
    LastSNP = ""
    min = 5  # Choix du nombre de lignée minimum données pour considérer le résultat intéressant
    total = 0
    for row in file_reader:
        if index == 0:
            headers = row
        elif row[0] != LastSNP:   # Enlever les SNP écrits 2 fois de suite avec gene ID !=
                                  # TO DO: Verifier que les lignes sont bien pareil??
            # Look for differences
            parse = dict()
            Nbrelignee = 0

            for i in range (8, len(row) - 1):
                if row[i] != '' and not "conflict" in row[i]:  # row[i].startswith("?"):
                                                               # Parfois ligne avec conflit sans forcément commmencer par un ?
                    Nbrelignee = Nbrelignee + 1
                    if row[i] in parse:
                        parse[row[i]].append(headers[i])
                    else:
                        parse[row[i]] = [ headers[i] ]
                # Check if there is a unique element - une lignée pour un allele donné
                # parse = dictionnaire des lignées pour chaque allele
                onlyOnce = True
                differenciator = ""
                if len(parse) > 1 and Nbrelignee > min:         # On élimine les cas où on a qu'un seul type d'allele
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
                result["SNP"] = row[0]
                result["diff"] = differenciator
                result["lignee"] = parse[differenciator][0]
                results.append(result)
                print (index, "SNP: ", row[0], " Differenciator: ", differenciator, " lignee: ", parse[differenciator][0], "Nombre de lignee: ", Nbrelignee)
                rst.write("\t".join([str(index), row[0], differenciator, parse[differenciator][0], str(Nbrelignee)]))
                rst.write("\n")
                total = total + 1
        LastSNP = row[0]
        index = index + 1

    # Use pandas library to build histogram
    df = pd.DataFrame(results)
    dt = pd.DataFrame(df.groupby('lignee').count())
    dt = dt.drop('diff', axis = 1)
    dt.columns = ['chr']
    print(dt)
    if limited_results.empty:
        limited_results = dt
    else:
        limited_results = pd.merge(left=limited_results, right=ft, on='lignee', how='outer')
    print("Total de SNP:    ", total)
    rst_output = open("Tableau recap chr lim.txt","w")
    rst_output.write(f"{dt}")
    rst_output.write("\n")
    rst_output.write("total de SNP" "\t" f"{total}" )
    rst_output.close()
    rst.close()
