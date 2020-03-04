import csv
import pandas as pd

with open("Mouse Chr lim.txt", "r") as file:
    rst = open("Result Chr lim.txt" , "w")
    rst.write("\t".join(["ligne", "gene", "allele", "lignee"]))
    rst.write("\n")
    file_reader = csv.reader(file, delimiter='\t')
    index = 0       # numéro de la ligne
    results =  []
    for row in file_reader:
        if index == 0:
            headers = row
        else:
            # CLook for differences
            parse = dict()
            for i in range (8, len(row) - 1):
                if row[i] != '' and not row[i].startswith("?"):
                    if row[i] in parse:
                        parse[row[i]].append(headers[i])
                    else:
                        parse[row[i]] = [ headers[i] ]
                # Check if there is a unique element - une lignée pour un allele donné
                # parse = dictionnaire des lignées pour chaque allele
                onlyOnce = True
                differenciator = ""
                if len(parse) > 1:          # On élimine les cas où on a qu'un seul type d'allele
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
                result["gene"] = row[0]
                result["diff"] = differenciator
                result["lignee"] = parse[differenciator][0]
                results.append(result)
                print (index, "gene: ", row[0], " Differenciator: ", differenciator, " lignee: ", parse[differenciator][0])
                rst.write("\t".join([str(index), row[0], differenciator, parse[differenciator][0]]))
                rst.write("\n")

        index = index + 1

    # Use pandas library to build histogram
    df = pd.DataFrame(results)
    dt = (df.groupby('lignee').count())
    print(dt)
    rst_output = open("Tableau recap chr lim.txt","w")
    rst_output.write( f"{dt}")
    rst_output.close()
    rst.close()
