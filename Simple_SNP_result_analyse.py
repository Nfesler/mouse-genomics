import csv
import pandas as pd

filename_pattern = "C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\Liste\\result Chr{0}.csv"
out_pattern = "C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\followed\\followed result Chr{0}.csv"

def analyse(filename, out):
    with open ("C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\Liste\\result Chr1.csv", 'r') as file:

        file_reader = csv.DictReader(file, delimiter = ',')
        index = 0
        last_strain = ''
        last_line = -2
        followed = False
        saved_line = -10
        results = []

        for row in file_reader:
            if index != 0:
                if int(row['line']) == (last_line + 1) and row['strain'] == last_strain:
                    followed = True
                    if int(row['line']) == (saved_line + 2):
                        # print(row)
                        results.append(row)
                    else:
                        # print(last_row, '\n', row)
                        results.append(last_row)
                        results.append(row)
                else:
                    followed = False

            if followed:
                saved_line = last_line

            last_row = row
            last_line = int(row['line'])
            last_strain = row['strain']
            index = index + 1

        dt = pd.DataFrame(results)
        dt = dt.drop( '' , axis = 1)
        dt.to_csv(out)

for i in range(19):
    index = i + 1
    filename = filename_pattern.format(index)
    out = out_pattern.format(index)
    print('Chr', index)
    analyse(filename, out)
