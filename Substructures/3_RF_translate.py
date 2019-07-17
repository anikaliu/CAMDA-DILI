import sys
import csv


def import_SMARTS():
    with open('/home/o/Desktop/smarts0.csv', newline='', mode='r') as csvfile:
        linereader = csv.reader(csvfile)  # iterable
        indices_to_fragments = {}
        for row in linereader:
            index = int(row[0])
            fragment = row[1]
            indices_to_fragments[index-1] = fragment
            # print(fragment)
    return indices_to_fragments


def import_RF_outcomes(indices_to_fragments):

    with open('/home/o/Desktop/indices0.csv', newline='', mode='r') as csvfile:
        linereader = csv.reader(csvfile)  # iterable
        for row in linereader:
            index = int(row[0])
            importance = row[1]
            fragment = indices_to_fragments[index]
            csv_row = [fragment, importance]
            joined = ','.join(csv_row)
            print(joined)


indices_to_fragments = import_SMARTS()
fragment, importance = import_RF_outcomes(indices_to_fragments)
