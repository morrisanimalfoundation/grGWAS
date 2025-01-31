import pandas as pd
import scipy.stats as stats
import sys

def main(file_path):
    # Load the headerless data
    data = pd.read_csv(file_path, sep='\t', header=None)

    # Assign column names
    data.columns = ['subject_id', 'binary_phenotype', 'gender']

    # Display the first few rows of the dataframe
    print("First few rows of the data:")
    print(data.head())

    # Create a contingency table
    contingency_table = pd.crosstab(data['binary_phenotype'], data['gender'])

    # Display the contingency table
    print("\nContingency Table:")
    print(contingency_table)

    # Perform the Chi-square test
    chi2, p, dof, ex = stats.chi2_contingency(contingency_table)

    # Display the results
    print(f"\nChi2 Statistic: {chi2}")
    print(f"P-value: {p}")

    if p < 0.05:
        print("There is a significant difference in gender status with respect to the binary phenotype.")
    else:
        print("There is no significant difference in gender status with respect to the binary phenotype.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <file_path>")
    else:
        file_path = sys.argv[1]
        main(file_path)

