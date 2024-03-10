import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import skew, kurtosis, normaltest

# Load the dataset
file_path = 'C:\\Users\\Musti Tanvir\\Desktop\\university work\\freelance work\\Projects\\Rafay\\Flood Prediction\\data\\kotri barrage (2).xlsx'
df = pd.read_excel(file_path, sheet_name='SUMMARY')  # Replace with your sheet name

# Function to perform comprehensive univariate analysis
def univariate_analysis(data, column):
    """
    Perform comprehensive univariate analysis for a given column.
    Includes descriptive statistics, histograms, boxplots, skewness, kurtosis, and normality test.
    """
    # Descriptive Statistics
    desc_stats = data[column].describe()

    # Skewness, Kurtosis, and Normality Test
    desc_stats['skewness'] = skew(data[column].dropna())
    desc_stats['kurtosis'] = kurtosis(data[column].dropna())
    desc_stats['normality_test_pval'] = normaltest(data[column].dropna())[1]

    # Visualization: Histogram and Boxplot
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    sns.histplot(data[column], bins=15, kde=True, ax=axes[0])
    axes[0].set_title(f'Histogram of {column}')

    sns.boxplot(y=data[column], ax=axes[1])
    axes[1].set_title(f'Boxplot of {column}')

    plt.tight_layout()
    plt.show()

    return desc_stats

# Perform univariate analysis for each numeric column and store results
analysis_results = {}
for column in df.select_dtypes(include=[np.number]).columns:
    analysis_results[column] = univariate_analysis(df, column)

# Convert results to a DataFrame in Excel
analysis_df = pd.DataFrame(analysis_results)

# Save the analysis results to an Excel file
output_file_path = 'kotri univariate.xlsx'
analysis_df.to_excel(output_file_path)
