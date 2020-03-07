import sys
import csv
import os
from collections import Counter

directory = sys.argv[1]

def extract_strain(file_name):
    with open(file_name, "r") as input_file:
        file_reader = csv.reader(input_file, delimiter = '\t')
        row = next(file_reader)
        strains = row[8:len(row) - 1]
        return strains

all_strains = []
for file_name in os.listdir(directory):
    if file_name.endswith(".txt"):
        all_strains.extend(extract_strain(file_name))

# Count strains occurrence
counter = Counter(all_strains)

# Find the maximum occurrence
occurencies = counter.values()
max_occurency = max(occurencies)

# filter all strain with max value
common_strains_dict = dict(filter(lambda el: el[1] >= max_occurency, counter.items()))

print(common_strains_dict.keys())
