import pandas as pd
from copulas.bivariate import Clayton, Frank, Gumbel
from scipy.stats import rankdata
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.stats as stats

# Load the dataset
file_path = 'path\\to\\your\\file' #Path to the Excel File
df = pd.read_excel(file_path, sheet_name='SUMMARY')

column1 = 'PD in m/s'
column2 = 'V in m/s'

# Transform data to uniform marginals
def to_uniform(data):
    ranks = rankdata(data, method='average')
    return (ranks - 1) / (len(data) - 1)

uniform_data1 = to_uniform(df[column1].dropna())
uniform_data2 = to_uniform(df[column2].dropna())

# Fit Bivariate Copulas
bivariate_clayton = Clayton()
bivariate_clayton.fit(np.column_stack([uniform_data1, uniform_data2]))

bivariate_frank = Frank()
bivariate_frank.fit(np.column_stack([uniform_data1, uniform_data2]))

bivariate_gumbel = Gumbel()
bivariate_gumbel.fit(np.column_stack([uniform_data1, uniform_data2]))

# Define the Joe Copula class
class Joe:
    def __init__(self, theta):
        self.theta = theta

    def pdf(self, u):
        n = u.shape[0]
        return np.prod((1 + (self.theta - 1) * u)**(-1/self.theta - 1), axis=1) * (1 - n)**(1/self.theta)

# Fit the Joe Copula
theta = 2  # Update with the desired value of theta
bivariate_joe = Joe(theta)

# Save the results to an Excel file (optional)
results_df = pd.DataFrame({
    f'{column1} (Uniform)': uniform_data1,
    f'{column2} (Uniform)': uniform_data2
})
results_df.to_excel('kotri bivariate.xlsx', index=False)

# Define the meshgrid for plotting
u = np.linspace(0.01, 0.99, 100)
v = np.linspace(0.01, 0.99, 100)
u, v = np.meshgrid(u, v)
uv = np.array([u.flatten(), v.flatten()]).T

# Plotting the results
# Clayton Copula
fig = plt.figure(figsize=(12, 8))

# 3D surface plot
ax = fig.add_subplot(231, projection='3d')
c = bivariate_clayton.pdf(uv).reshape((100, 100))
ax.plot_surface(u, v, c, cmap='viridis')
ax.set_title('Clayton Copula 3D Surface Plot')

# Scatter plot
ax = fig.add_subplot(232)
ax.scatter(uniform_data1, uniform_data2)
ax.set_title('Clayton Copula Scatter Plot')

# QQ plot
ax = fig.add_subplot(233)
stats.probplot(uniform_data1, dist="uniform", plot=ax)
ax.set_title('Clayton Copula QQ Plot')

# Frank Copula
# 3D surface plot
ax = fig.add_subplot(234, projection='3d')
c = bivariate_frank.pdf(uv).reshape((100, 100))
ax.plot_surface(u, v, c, cmap='viridis')
ax.set_title('Frank Copula 3D Surface Plot')

# Scatter plot
ax = fig.add_subplot(235)
ax.scatter(uniform_data1, uniform_data2)
ax.set_title('Frank Copula Scatter Plot')

# QQ plot
ax = fig.add_subplot(236)
stats.probplot(uniform_data1, dist="uniform", plot=ax)
ax.set_title('Frank Copula QQ Plot')

plt.tight_layout()
plt.show()

# Gumbel Copula
fig = plt.figure(figsize=(12, 8))

# 3D surface plot
ax = fig.add_subplot(231, projection='3d')
c = bivariate_gumbel.pdf(uv).reshape((100, 100))
ax.plot_surface(u, v, c, cmap='viridis')
ax.set_title('Gumbel Copula 3D Surface Plot')

# Scatter plot
ax = fig.add_subplot(232)
ax.scatter(uniform_data1, uniform_data2)
ax.set_title('Gumbel Copula Scatter Plot')

# QQ plot
ax = fig.add_subplot(233)
stats.probplot(uniform_data1, dist="uniform", plot=ax)
ax.set_title('Gumbel Copula QQ Plot')

plt.tight_layout()
plt.show()

# Joe Copula
fig = plt.figure(figsize=(12, 8))

# 3D surface plot
ax = fig.add_subplot(231, projection='3d')
c = bivariate_joe.pdf(uv).reshape((100, 100))
ax.plot_surface(u, v, c, cmap='viridis')
ax.set_title('Joe Copula 3D Surface Plot')

# Scatter plot
ax = fig.add_subplot(232)
ax.scatter(uniform_data1, uniform_data2)
ax.set_title('Joe Copula Scatter Plot')

# QQ plot
ax = fig.add_subplot(233)
stats.probplot(uniform_data1, dist="uniform", plot=ax)
ax.set_title('Joe Copula QQ Plot')

plt.tight_layout()
plt.show()
