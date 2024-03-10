import pandas as pd
from copulas.univariate import GaussianKDE
from copulas.bivariate import Clayton
from scipy.stats import rankdata
import numpy as np
import plotly.graph_objects as go


# Load the dataset
file_path = 'C:\\Users\\Musti Tanvir\\Desktop\\university work\\freelance work\\Projects\\Rafay\\Flood Prediction\\data\\kotri barrage (2).xlsx'
df = pd.read_excel(file_path, sheet_name='SUMMARY')

column1 = 'PD in m/s'  # Update with actual column name
column2 = 'V in m/s'   # Update with actual column name

# Transform data to uniform marginals
def to_uniform(data):
    ranks = rankdata(data, method='average')
    return (ranks - 1) / (len(data) - 1)

uniform_data1 = to_uniform(df[column1].dropna())
uniform_data2 = to_uniform(df[column2].dropna())

# Fit Univariate and Bivariate Copulas
univariate_copula1 = GaussianKDE()
univariate_copula1.fit(uniform_data1)

univariate_copula2 = GaussianKDE()
univariate_copula2.fit(uniform_data2)

bivariate_copula = Clayton()
bivariate_copula.fit(np.column_stack([uniform_data1, uniform_data2]))

# Save the results to an Excel file
results_df = pd.DataFrame({
    f'{column1} (Uniform)': uniform_data1,
    f'{column2} (Uniform)': uniform_data2
})
results_df.to_excel('kotri bivariate.xlsx', index=False)


# Generate samples for visualization
sampled_data = bivariate_copula.sample(500)

# Plotting the results with Plotly
fig = go.Figure()
fig.add_trace(go.Scatter(x=uniform_data1, y=uniform_data2, mode='markers', name='Original Data'))
fig.add_trace(go.Scatter(x=sampled_data[:, 0], y=sampled_data[:, 1], mode='markers', name='Sampled Data'))
fig.update_layout(title='Bivariate Copula Analysis', xaxis_title=column1, yaxis_title=column2)

# Display or save the plot as per your requirement
fig.show()
