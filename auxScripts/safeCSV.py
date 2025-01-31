import csv
import sys

def process_csv(input_file):
    with open(input_file, 'r', newline='') as infile:
        reader = csv.reader(infile, quotechar='"', delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer = csv.writer(sys.stdout, delimiter=',', quoting=csv.QUOTE_NONE, escapechar='\\', lineterminator='\n')

        for row in reader:
            new_row = [field.replace(',', ';').replace('\n', ';') for field in row]
            writer.writerow(new_row)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    process_csv(input_file)


