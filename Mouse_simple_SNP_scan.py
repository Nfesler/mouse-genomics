import csv
import pandas as pd

filename_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Génome souris\\Mouse Chr {0}.txt'
out1_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\Liste\\result Chr{0}.txt'
out2_pattern = 'C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\StrainTab\\result 2 Chr{0}.txt'

full_results = pd.DataFrame()

def analyse(filename, out1, out2, columnName):
    with open(filename, "r") as file:
        rst = open(out1 , "w")
        rst.write("\t".join(["index", "SNP", "allele", "lignee", "Nbre de lignee"]))
        rst.write("\n")
        file_reader = csv.reader(file, delimiter='\t')
        index = 0       # numéro de la ligne
        results =  []
        LastSNP  = ""
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
                    # print (index, "SNP: ", row[0], " Differenciator: ", differenciator, " lignee: ", parse[differenciator][0], "Nombre de lignee: ", Nbrelignee)
                    rst.write("\t".join([str(index), row[0], differenciator, parse[differenciator][0], str(Nbrelignee)]))
                    rst.write("\n")
                    total = total + 1
                index = index + 1
            LastSNP = row[0]

        # Use pandas library to build histogram
        global full_results
        df = pd.DataFrame(results)
        dt = pd.DataFrame(df.groupby('lignee').count())
        dt = dt.drop('diff', axis=1)
        dt.columns = [columnName]
        print(dt)
        print("Total de SNP:    ", total)
        if full_results.empty:
            full_results = dt
        else:
            full_results = pd.merge(left=full_results, right=dt, on='lignee', how='outer')
        rst_output = open(out2,"w")
        rst_output.write(f"{dt}\n")
        rst_output.write("Total de SNPs:" "\t" f"{total}")
        rst_output.close()
        rst.close()

for i in range(19):
    index = i + 1
    filename = filename_pattern.format(index)
    out1 = out1_pattern.format(index)
    out2 = out2_pattern.format(index)
    print("chromosome", index)
    analyse(filename, out1, out2, columnName='chr' + str(index))

print(full_results)
full_results.to_csv('C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\full_results_simple.csv')
