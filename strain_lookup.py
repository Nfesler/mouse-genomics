import sys
import csv
import os
from collections import Counter

#directory = sys.argv[1]
directory = 'Génome souris'

def extract_strain(file_name):
    with open(file_name, "r") as input_file:
        file_reader = csv.reader(input_file, delimiter = '\t')
        row = next(file_reader)
        strains = row[8:len(row) - 1]
        return strains

all_strains = []
for file in os.scandir(directory):
    if file.path.endswith(".txt"):
        all_strains.extend(extract_strain(file.path))

# Count strains occurrence
counter = Counter(all_strains)
# print(counter)

# Find the maximum occurrence
occurencies = counter.values()
max_occurency = max(occurencies)
#print(max_occurency)

# filter all strain with max value
common_strains_dict = dict(filter(lambda el: el[1] >= max_occurency, counter.items()))
# print(common_strains_dict)

print(common_strains_dict.keys())
