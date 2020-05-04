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
    rst = open("Result Windowed Chr lim.txt" , "w")
    file_reader = csv.reader(file, delimiter='\t')
    index = 0       # numéro de la ligne
    LastSNP = ""
    min_strain = 5  # Choix du nombre de lignée minimum données pour considérer le résultat intéressant
    total = 0
    results = pd.DataFrame(columns=['line', 'snp', 'strain', 'n_samples'])

    for row in file_reader:
        if index == 0:
            headers = row
            index = index + 1
        elif row[0] != LastSNP:   # Enlever les SNP écrits 2 fois de suite avec gene ID !=

            line_index = (index - 1 ) % WINDOW_SIZE

            strain_data[line_index] = np.array((row[8:len(row) - 1])) # La dernière colonne n'a pas une valeur
            snp_pos[line_index] = row[0]

            # print('\n')
            # print(index)
            # print(strain_data)


            if index > (WINDOW_SIZE - 1) :
                stacked_strain_data = np.stack(strain_data)
                # Start calculation
                dbscan_data = build_dbscan_data(stacked_strain_data)
                #print(dbscan_data)
                dbscan_result = DBSCAN(metric=strain_comparator,eps=5,min_samples=2,algorithm='brute').fit(dbscan_data)
                labels = dbscan_result.labels_

                n_unique = np.where(labels < 0)

                # print(labels)
                # print(n_unique)
                if (len(n_unique[0]) == 1):
                    data_index = n_unique[0][0]
                    # print(data_index)
                    strain_index = dbscan_data[data_index][0]
                    # print(strain_index)
                    # print(index - 1)
                    # print(snp_pos)
                    # print(snp_pos[(index) % WINDOW_SIZE])
                    # print(headers[strain_index + 8])
                    # print(len(labels))
                    results = results.append({'line': index -1, 'snp': snp_pos[(index) % WINDOW_SIZE], 'strain':headers[strain_index + 8], 'n_samples': len(labels)}, ignore_index = True)
                    total = total + 1

            index = index + 1

        LastSNP = row[0]

    print(results)
    rst.write(str(results))
    rst.close()
    dt = pd.DataFrame(results.groupby('strain').count()[['snp']])
    print(dt)
    print('Total:', total)
