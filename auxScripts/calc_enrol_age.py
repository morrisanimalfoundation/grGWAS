import csv
import sys
from datetime import datetime

def process_csv(input_file):
    with open(input_file, 'r', newline='') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['enrollment_age']
        
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter=',', quoting=csv.QUOTE_NONE, escapechar='\\', lineterminator='\n')
        writer.writeheader()

        for row in reader:
            try:
                birth_date = datetime.strptime(row['birth_date'], '%Y-%m')
                last_reported_date = datetime.strptime(row['enrolled_date'], '%Y-%m')
                difference_in_days = (last_reported_date - birth_date).days
                difference_in_years = round(difference_in_days / 365, 2)
            except ValueError as e:
                difference_in_years = 'N/A'  # In case of any date parsing errors

            row['enrollment_age'] = difference_in_years
            writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    process_csv(input_file)

