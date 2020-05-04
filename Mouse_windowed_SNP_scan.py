import csv
import sys
import pandas as pd
from sklearn.cluster import DBSCAN
import numpy as np
from progress.bar import ShadyBar

filename_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Génome souris\\Mouse Chr {0}.txt'
out_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Windowed\\Liste\\result windowed chr {0}.csv'

full_results = pd.DataFrame()

WINDOW_SIZE = 3

MIN_SAMPLES = 3

MAXIMUM_DISTANCE = 1000000

strain_data = np.zeros(WINDOW_SIZE, dtype = 'object')

SNP_pos = np.zeros(WINDOW_SIZE, dtype = 'object')





def analyse(filename, out, nb_chr):

    def line_count(filename):
        with open (filename) as file:
            file_reader = csv.reader(file)
            for i, l in enumerate(file_reader):
                pass
            return (i + 1)

    def simple_comparator(x,y):
        # Comparaison des SNP de 2 lignée un par rapport à l'autre
        if x =='' or y == '' :
            # empty SNP could be anything
            return 0
        elif x == y :
            # Si 2 SNP sont identiques
            return 0
        else :
            return MAXIMUM_DISTANCE


    def strain_comparator(x,y):
        # comparaison de l'ensemble de la fenêtre
        i, j = int(x[0]), int(y[0])
        distance = 0
        sdt_i = stacked_strain_data[:,i]
        sdt_j = stacked_strain_data[:,j]

        if (sdt_i == '').all() or (sdt_j == '').all():
            # Enlever si les 2 lignée sont vides (deja fait lors de simple_validator)
            distance = MAXIMUM_DISTANCE
        else :
            for k in range(WINDOW_SIZE):
                distance = min(distance + simple_comparator(stacked_strain_data[k][i], stacked_strain_data[k][j]), MAXIMUM_DISTANCE)

        return distance


    def simple_validator(data):
        search_conflict = np.core.defchararray.find(data, 'conflict') != -1
            # si conflict dans une valeur alors valeur = -1 donc valeur = False
        if search_conflict.any():
            # Si présence d'un False dans le tableau
            return False
        elif (data == '').all():
            # Si tout les SNP sont vides
            return False
        else:
            return True


    def build_dbscan_data(st_data):
        dbscan_data = []
        for i in range(len(st_data[0])):
            if simple_validator(st_data[:,i]):
                dbscan_data = dbscan_data + [[i]]
        return dbscan_data


    with open(filename, 'r') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        total_lines = line_count(filename)
        index = 0
        lastSNP = ''
        total = 0
        results = pd.DataFrame(columns = ['Line', 'SNP', 'Strain', 'n_samples'])
        global full_results
        bar = ShadyBar(nb_chr, max=total_lines / 1000, suffix='%(percent)d%%')


        for row in file_reader:
            if index == 0:
                headers = row
                index = index + 1
            elif row[0] != lastSNP:
                # Enlever les SNP écrits 2 fois de suite avec gene ID !=
                line_index = (index - 1) % WINDOW_SIZE
                # Reste de la division de index par windows_size
                strain_data[line_index] = np.array((row[8:len(row) - 1])) # la dernière colonne n'a pas de valeur
                # Sur le tableau strain_data, la ligne line_index est remplacée par la ligne correspondante
                SNP_pos[line_index] = row[0]
                # Enregistrer le SNP

                if index > (WINDOW_SIZE - 1) : # Eviter les 2 première lignes
                    stacked_strain_data = np.stack(strain_data) # creer tableau a 2D a partir de tableau de ligne
                    # Start calculation:
                    dbscan_data = build_dbscan_data(stacked_strain_data)
                    if len(dbscan_data) >= MIN_SAMPLES:
                        dbscan_result = DBSCAN(metric = strain_comparator, eps = 5, min_samples = 2, algorithm = 'brute').fit(dbscan_data)
                        # Assembler les strain en groupe identiques selon la technique de strain_comparator
                        # Groupe 0; 1; 2; ... et -1 si seul
                        labels = dbscan_result.labels_
                        n_unique = np.where(labels < 0)
                        # Retrouver les pattern unique qui ont une valeur de -1

                        if(len(n_unique[0]) == 1):
                            data_index = n_unique[0][0]
                            strain_index = dbscan_data[data_index][0]
                            results = results.append({'Line': index - 1, 'SNP': SNP_pos[(index) % WINDOW_SIZE], 'Strain': headers[strain_index + 8], 'n_samples': len(labels)}, ignore_index = True)
                            total = total + 1


                if index % 1000 == 0:
                    bar.next()
                    # progress = round((index / total_lines) * 100, 2)
                    # print(progress, '% ', end='\r')

                index = index + 1

            lastSNP = row[0]

        bar.finish()

        print('\n')
        print(results)
        rst.to_csv(out)
        dt = pd.DataFrame(results.groupby('Strain').count()[['SNP']])
        print(dt)
        print('Total de SNP:    ', total)
        dt.columns = [nb_chr]
        if full_results.empty :
            full_results = dt
        else :
            full_results = pd.merge(left = full_results, right = dt, on = 'Strain', how = 'outer')



for i in range(19):
    index = i + 1
    filename = filename_pattern.format(index)
    out = out_pattern.format(index)
    # print("chromosome", index)
    analyse(filename, out, nb_chr = 'chr' + str(index))

print(full_results)
full_results.to_csv('C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Windowed\\full_results_windowed.csv')
