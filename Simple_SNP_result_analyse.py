import csv

filename_pattern = "C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\Liste\\result Chr{0}.txt"

def analyse(filename):
    with open ("C:\\Users\\nicol\\OneDrive\\Documents\\GitHub\\mouse-genomics\\Result\\Simple\\Liste\\result Chr1.txt", 'r') as file:
        file_reader = csv.reader(file, delimiter = '\t')
        index = 0
        line = -2
        for row in file_reader:
            if index == 0:
                headers = row
            else:
                if int(row[0]) == (line + 1):
                    print(line, 'true')
                line = int(row[0])
            index = index + 1

for i in range(19):
    index = i + 1
    filename = filename_pattern.format(index)
    print('Chr', index)
    analyse(filename)
