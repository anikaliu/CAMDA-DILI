import sys
import csv

def import_dataset():
    with open('fragments_and_compounds_input.csv', newline='', mode='r') as csvfile:
        linereader = csv.reader(csvfile)  # iterable
        fragments_to_compounds = {}
        for row in linereader:
            fragment = row[0]
            compounds = row[1:]
            compounds = [compound for compound in compounds if compound != '?']
            fragments_to_compounds[fragment] = compounds
            # print(fragment)
    return(fragments_to_compounds)

def import_datalabels():
    with open('compounds_MOSS input.csv', newline='', mode='r') as csvfile:
        linereader = csv.reader(csvfile)  # iterable
        compound_to_label = {}
        next(linereader)  # skip the file header - the 1st row
        for row in linereader:
            compound = row[0]
            label = row[2]
            compound_to_label[compound] = label
            # print(fragment)
    return compound_to_label

fragments_to_compounds = import_dataset()
compounds_to_fragments = {}
compound_to_label = import_datalabels()
all_fragments = []
for fragment, compounds in fragments_to_compounds.items():
    for compound in compounds:
        if not compound in compounds_to_fragments:
            compounds_to_fragments[compound] = []
        compounds_to_fragments[compound].append(fragment)
        if not fragment in all_fragments:
            all_fragments.append(fragment)

header = ['DILI Label', 'Compound Name', 'Fragment Count'] + all_fragments
print(','.join(header))
for compound, fragments in compounds_to_fragments.items():
    label = compound_to_label[compound]
    binary_vector = []

    fragment_count = 0
    for fragment in all_fragments:

        if fragment in fragments:
            binary_vector.append('1')
            fragment_count += 1
        else:
            binary_vector.append('0')
    csv_row = [label, compound, str(fragment_count)]+binary_vector
    joined = ','.join(csv_row)
    print(joined)
