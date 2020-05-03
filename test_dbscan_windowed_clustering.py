import csv
import sys
import pandas as pd
from sklearn.cluster import DBSCAN
import numpy as np

limited_results = pd.DataFrame()

WINDOW_SIZE = 3

MAXIMUM_DISTANCE = 1000000

strain_data = np.zeros(WINDOW_SIZE, dtype= 'object')

snp_pos = np.zeros(WINDOW_SIZE, dtype= 'object')



def simple_comparator(x,y):
    if x == '' or y == '':
        # empty could be anything
        return 0
    elif x == y:
        return 0
    else :
        return MAXIMUM_DISTANCE

def strain_comparator(x, y):
    i, j = int(x[0]), int(y[0])

    distance = 0

    # If none information is present, return large distance
    sdt_i = stacked_strain_data[:,i]
    sdt_j = stacked_strain_data[:,j]

    if (sdt_i == '').all() or (sdt_j == '').all():
        distance = MAXIMUM_DISTANCE
    else :
        for k in range(WINDOW_SIZE):
            distance = min(distance + simple_comparator(stacked_strain_data[k][i], stacked_strain_data[k][j]), MAXIMUM_DISTANCE)

    return distance

def simple_validator(data):
    search_conflict = np.core.defchararray.find(data, 'conflict') != -1
    if search_conflict.any():
        return False
    elif (data == '').all():
        return False
    else:
        return True

def build_dbscan_data(st_data):
    # We eliminate any strain with 'conflict' or a complete blank window
    dbscan_data = []
    for i in range(len(st_data[0])):
        if simple_validator(st_data[:,i]):
            dbscan_data = dbscan_data + [[i]]
    return dbscan_data




with open("Mouse Chr lim.txt", "r") as file:
    rst = open("Result Chr lim.txt" , "w")
    rst.write("\t".join(["ligne", "SNP", "allele", "lignee", "Nbre de lignee"]))
    rst.write("\n")
    file_reader = csv.reader(file, delimiter='\t')
    index = 0       # numéro de la ligne
    LastSNP = ""
    min_strain = 5  # Choix du nombre de lignée minimum données pour considérer le résultat intéressant
    total = 0
    results = pd.DataFrame(columns=['line', 'snp', 'strain'])

    for row in file_reader:
        if index == 0:
            headers = row
            index = index + 1
        elif row[0] != LastSNP:   # Enlever les SNP écrits 2 fois de suite avec gene ID !=

            line_index = (index - 1 ) % WINDOW_SIZE

            strain_data[line_index] = np.array((row[8:len(row) - 1])) # La dernière colonne n'a pas une valeur
            snp_pos[line_index] = row[0]

            #print('\n')
            #print(index)
            print(strain_data)


            if index > (WINDOW_SIZE - 1) :
                stacked_strain_data = np.stack(strain_data)
                # Start calculation
                dbscan_data = build_dbscan_data(stacked_strain_data)
                #print(dbscan_data)
                dbscan_result = DBSCAN(metric=strain_comparator,eps=5,min_samples=2,algorithm='brute').fit(dbscan_data)
                labels = dbscan_result.labels_

                n_unique = np.where(labels < 0)

                #print(labels)
                if (len(n_unique[0]) == 1):
                    data_index = n_unique[0][0]
                    strain_index = dbscan_data[data_index][0]
                    print(index - 1)
                    print(snp_pos)
                    print(snp_pos[(index) % WINDOW_SIZE])
                    print(headers[strain_index + 8])
                    results = results.append({'line': index -1, 'snp': snp_pos[(index) % WINDOW_SIZE], 'strain':headers[strain_index + 8]}, ignore_index=True)

            #
            #
            # # Look for differences
            # parse = dict()
            # Nbrelignee = 0
            # for i in range (8, len(row) - 1):
            #     if row[i] != '' and not "conflict" in row[i]:  # row[i].startswith("?"):
            #                                                    # Parfois ligne avec conflit sans forcément commmencer par un ?
            #         Nbrelignee = Nbrelignee + 1
            #         if row[i] in parse:
            #             parse[row[i]].append(headers[i])
            #         else:
            #             parse[row[i]] = [ headers[i] ]
            #     # Check if there is a unique element - une lignée pour un allele donné
            #     # parse = dictionnaire des lignées pour chaque allele
            #     onlyOnce = True
            #     differenciator = ""
            #     if len(parse) > 1 and Nbrelignee > min:         # On élimine les cas où on a qu'un seul type d'allele
            #                                                     # ou un nombre insuffisant de lignées pour que ce soit
            #                                                     # intéressant
            #         for k,v in parse.items():
            #             if len(v) == 1:
            #                 if differenciator == "":
            #                     differenciator = k
            #                 else:
            #                     # La deuxième fois qu'on rencontre un allele avec une seule
            #                     # lignée, differenciator <> "" donc on positionne le flag onlyOnce à False
            #                     # car cela veut dire qu'on a plusieurs alleles avec une seule lignée
            #                     onlyOnce = False
            #                     break
            #
            # if onlyOnce and (differenciator != ""):
            #     result = dict()
            #     result["SNP"] = row[0]
            #     result["diff"] = differenciator
            #     result["lignee"] = parse[differenciator][0]
            #     results.append(result)
            #     print (index, "SNP: ", row[0], " Differenciator: ", differenciator, " lignee: ", parse[differenciator][0], "Nombre de lignee: ", Nbrelignee)
            #     rst.write("\t".join([str(index), row[0], differenciator, parse[differenciator][0], str(Nbrelignee)]))
            #     rst.write("\n")
            #     total = total + 1
            index = index + 1
        LastSNP = row[0]


    # Use pandas library to build histogram
    # df = pd.DataFrame(results)
    # dt = pd.DataFrame(df.groupby('lignee').count())
    # dt = dt.drop('diff', axis = 1)
    # dt.columns = ['chr']
    # print(dt)
    # if limited_results.empty:
    #     limited_results = dt
    # else:
    #     limited_results = pd.merge(left=limited_results, right=ft, on='lignee', how='outer')
    # print("Total de SNP:    ", total)
    # rst_output = open("Tableau recap chr lim.txt","w")
    # rst_output.write(f"{dt}")
    # rst_output.write("\n")
    # rst_output.write("total de SNP" "\t" f"{total}" )
    # rst_output.close()
    rst.close()
    print(results)
    dt = pd.DataFrame(results.groupby('strain').count()[['snp']])
    print(dt)
