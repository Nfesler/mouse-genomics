import csv
import pandas as pd

with open("Mouse Chr lim.txt", "r") as file:
    file_reader = csv.reader(file, delimiter='\t')
    index = 0
    results =  []
    for row in file_reader:
        if index == 0:
            headers = row
        else:
            # CLook for differences
            parse = dict()
            for i in range (8, len(row) - 1):
                if row[i] <> '' and not row[i].startswith("?"):
                    if row[i] in parse:
                        parse[row[i]].append(headers[i])
                    else:
                        parse[row[i]] = [ headers[i] ]
                # Check if there is a unique element
                onlyOnce = True
                differenciator = ""
                if len(parse) > 1:
                    for k,v in parse.items():
                        if len(v) == 1:
                            if differenciator == "":
                                differenciator = k
                            else:
                                onlyOnce = False

            if onlyOnce and (differenciator != ""):
                result = dict()
                result["gene"] = row[0]
                result["diff"] = differenciator
                result["lignee"] = parse[differenciator][0]
                results.append(result)
                # print "gene: ", row[0], " Differenciator: ", differenciator, " lignee: ", parse[differenciator][0]

        index = index + 1
        
    df = pd.DataFrame(results)
    print(df.groupby('lignee').count())
