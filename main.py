# Libraries
import numpy as np
# from scipy.stats import norm
from tqdm import tqdm
import scipy
import scipy.stats
import csv
import subprocess
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.animation import FuncAnimation
from matplotlib.animation import ArtistAnimation
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image
from IPython.display import HTML

# Functions
from functions import get_csv_line, data, sum_mutants_gene, normalization, gaussian_distribution, proportion_in_interval, gaussian_mixture, norm_cdf, y_prop, metropolis_hastings, simulation_mutants, sampling_mutants, simulated_cml_hist

##################################### Data #################################

# Number of processed row for a given treatment (the number of tested genes)
n_row = 11623
# n_row = 50

# Path to files
treatment_name = 'data_example'
csv_file = f'{treatment_name}.csv'
print(f'\n{treatment_name}.csv processing ...')
table_list = data(csv_file,n_row)

# Histogram data - Sum lines for 2 x 9 contingency tables
data_example_sum_line_2x9 = [84, 253, 50, 1152, 33463, 1340, 748, 358, 97]

# Linear values for Phi
x_values_linear = [1e-9, 0.062, 0.12, 0.25, 0.5, 2, 4, 8, 16, 32]
x_values_center = [0.062/2, (0.062+0.12)/2, (0.12+0.25)/2, (0.25+0.5)/2, (0.5+2)/2, 3, 6, 12, 24]
x_values_continuous = np.linspace(x_values_linear[0],x_values_linear[len(x_values_linear)-1], 1000)

# Log values for Phi
x_values = np.log(x_values_linear)
x_values_center_log = np.log(x_values_center)
x_values_continuous_log = np.log(x_values_continuous)

# Choice of the studied treatment
treatment = data_example_sum_line_2x9
normalized_treatment = normalization(treatment)

# Target histogram data and center values
target_hist = np.array(normalized_treatment)
center_values = np.array(x_values_center_log)

################### Settings for MH algorithm estimation ###################

# Parameters for a given treatment
mu_0 = 0 # To fix for each treatment
sigma_0 = 0.01 # To fix for each treatment
print("\nmu_0:", mu_0)
print("sigma_0:", sigma_0)
print("\n")
num_samples = 10000 # Increased the number of samples for better convergence

# Gibbs sampler
gibbs = 1 # activation of Gibbs sampler in MH algorithm
"""if gibbs == 0:
    print('MH run without Gibbs sampler')
else:
    print('MH run with Gibbs sampler')"""

# Gap to the gaussian mixture model parameter sigma_m (if not activated it is equal to 1)
var_m = 1e-8 # Gibbs sampler activated
# var_m = 1e-2 # Gibbs sampler not activated
sigma_m = np.sqrt(var_m)

# Initialization
initial_p = proportion_in_interval(target_hist, center_values)
initial_parameters = [initial_p, np.mean(center_values), np.std(center_values)]

# Fixed or not fixed parameters without Gibbs sampler
# fixed_p = 'None'
fixed_p = 'None'
fixed_mu_1 = 'None'
fixed_sigma_1 = 'None'

# Run MH algorithm
samples = metropolis_hastings(normalized_treatment, x_values, num_samples, mu_0, sigma_0, initial_parameters,var_m, gibbs, fixed_p, fixed_mu_1, fixed_sigma_1)

## Extract the estimated parameters from MH algorithm
estimated_p, estimated_mu_1, estimated_sigma_1 = np.mean(samples, axis=0)
print("sigma_m setting: ", sigma_m)
print("estimated p: ", estimated_p)
print("estimated mu_1: ", estimated_mu_1)
print("estimated sigma_1: ", estimated_sigma_1)

############################### Simulations ###############################

# Simulate mutants according to estimated gaussian mixture parameters
z_sampling, simulated_table_list = simulation_mutants(table_list, n_row, mu_0, sigma_0, estimated_mu_1, estimated_sigma_1, estimated_p)

# Convert log values in linear values
simulated_table_list_exp = [[np.exp(value) for value in sublist] for sublist in simulated_table_list]

# Final simulation with round data
simulated_treatment_data = sampling_mutants(x_values, simulated_table_list)

# Specify the output CSV file path
output_csv_file = f'{treatment_name}_simulated.csv'

# Read the existing CSV file and extract the first two columns
existing_data = []
with open(csv_file, 'r', newline='') as file:
    reader = csv.reader(file, delimiter='\t')  # Use tab as the delimiter
    header = next(reader)
    for row in reader:
        existing_data.append(row[:2])

# Define the header for the new CSV file
header = ['gene']+[f'alleles_with_growth_{x}-{y}' for x, y in zip(x_values_linear[:-1], x_values_linear[1:])]

# Combine the header and data to form rows for the new CSV
csv_data = [header]
for i in range(len(simulated_treatment_data)):
    # Only include the simulated treatment data for each row
    row = existing_data[i]
    row = row[0].split(';')
    row = [int(x) if x.isdigit() else x for x in row]
    row = row[:1]  # Sélectionne les deux premiers éléments
    row = row + simulated_treatment_data[i]
    csv_data.append(row)

# Write the data to the new CSV file
with open(output_csv_file, 'w', newline='') as file:
    writer = csv.writer(file, delimiter=',')  # Use comma as the delimiter
    writer.writerow(header)
    writer.writerows(csv_data[1:])

print(f'\nCSV file "{output_csv_file}" has been created.')

# Run R code to perform Fisher's exact code
subprocess.run(['Rscript', 'Fisher-exact-test.R'])
print(f'\nCSV file "{treatment_name}_simulated_p_values.csv" has been created by running "Fisher-exact-test.R".')

# Cumulative histogram of simulated data
cml_hist = simulated_cml_hist(n_row, simulated_treatment_data)

############################ Results and graph ############################

## Treatment histogram and simulated histogram (pb continuous gaussian)
print("\ntreatment:", treatment)
print("simulated treatment:", cml_hist)
print("\n")

estimated_mixture = gaussian_mixture(x_values_continuous_log, estimated_p, mu_0, sigma_0, estimated_mu_1, estimated_sigma_1) / np.sum(gaussian_mixture(x_values_continuous_log, estimated_p, mu_0, sigma_0, estimated_mu_1, estimated_sigma_1))
# estimated_gaussian_0 = estimated_p * gaussian_distribution(x_values_continuous_log, mu_0, sigma_0)
# estimated_gaussian_1 = (1 - estimated_p) * gaussian_distribution(x_values_continuous_log, estimated_mu_1, estimated_sigma_1)

# Calcul du facteur de multiplication pour obtenir estimated_mixture non normalisé
total_count_cml_hist = sum(cml_hist)  # cml_hist is not normalized
estimated_mixture_non_normalized_cml_hist = estimated_mixture * total_count_cml_hist

"""
plt.plot(x_values_continuous_log, estimated_mixture_non_normalized_cml_hist, color = 'blue', label='Gaussian mixture model')
plt.bar(x_values_center_log, treatment, alpha=0.5, color = 'red', label=f'{treatment_name} experimental data')
plt.bar(x_values_center_log, cml_hist, alpha=0.5, color = 'blue', label=f'{treatment_name} simulated data')
plt.legend()
plt.xlabel('$\log{\Phi}$')
plt.ylabel('mutants count')
plt.title('Treatment and simulated Histogram')
plt.show()"""

# Test on non normalized histogram
estimated_mixture = gaussian_mixture(x_values_continuous_log, estimated_p, mu_0, sigma_0, estimated_mu_1, estimated_sigma_1)
plt.plot(x_values_continuous_log, estimated_mixture, color = 'blue', label='Gaussian mixture model')
plt.bar(x_values_center_log, treatment/np.sum(treatment), alpha=0.5, color = 'red', label=f'{treatment_name} experimental data')
plt.bar(x_values_center_log, cml_hist/np.sum(cml_hist), alpha=0.5, color = 'blue', label=f'{treatment_name} simulated data')
plt.legend()
plt.xlabel('$\log{\Phi}$')
plt.ylabel('mutants count')
plt.title('Treatment and simulated Histogram')
plt.show()

# Test for log2 in x values
# x_values_log2 = [2 ** val for val in x_values_center]
# plt.xscale('log')

## Parameters histogram
liste_p = samples[:,0]
liste_mu_1 = samples[:,1]
liste_sigma_1 = samples[:,2]

plt.hist(liste_p, range = (0, 1), bins = 20, color = 'red', alpha = 0.7, label ='$p$')
plt.hist(liste_mu_1, range = (-np.log(8), np.log(8)), bins = 20, color = 'blue', alpha = 0.7, label ='$\mu_1$')
plt.hist(liste_sigma_1, range = (0, np.log(16)), bins = 20, color = 'orange', alpha = 0.7, label ='$\sigma_1$')
plt.title('Histograms of gaussian mixture parameters')
plt.xlabel('$\log{\Phi}$')
plt.ylabel('count')
plt.grid()
plt.legend()
plt.show()

## Parameters according to iterations
iterations = np.arange(0,num_samples,1)

plt.subplot(3, 1, 1)
plt.plot(iterations, liste_p, color = 'red', alpha = 0.7, label =f'$E[p] = {estimated_p:.2f}$', linewidth=0.8)
plt.title('Evolution of $p$')
plt.xlabel('iterations')
plt.ylabel('parameter value')
plt.grid()
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(iterations, liste_mu_1, color = 'blue', alpha = 0.7, label =f'$E[\mu_1] = {estimated_mu_1:.2f}$', linewidth=0.8)
plt.title('Evolution of $\mu_1$')
plt.xlabel('iterations')
plt.ylabel('parameter value')
plt.grid()
plt.legend()

plt.subplot(3, 1, 3)
plt.plot(iterations, liste_sigma_1, color = 'orange', alpha = 0.7, label =f'$E[\\sigma_1] = {estimated_sigma_1:.2f}$', linewidth=0.8)
plt.title('Evolution of $\sigma_1$')
plt.xlabel('iterations')
plt.ylabel('parameter value')
plt.grid()
plt.legend()

plt.tight_layout()
plt.subplots_adjust(hspace=0.7)
plt.show()

## Parameter space in 3D graph
parametres = samples
iterations = np.arange(num_samples)

# Séparation des paramètres en trois listes pour chaque dimension
liste_p = parametres[:, 0]
liste_mu_1 = parametres[:, 1]
liste_sigma_1 = parametres[:, 2]

# Création de la figure 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Représentation des points en 3D
ax.scatter(liste_p, liste_mu_1, liste_sigma_1, c='black', marker='o', s=5, alpha = 0.7)

# Dégradé de couleur selon les itérations (utilisation d'une carte de couleur "plasma" pour illustrer)
sc = ax.scatter(liste_p, liste_mu_1, liste_sigma_1, c=iterations, cmap='plasma', marker='o', s=10, alpha=0.7)

## TODO: Add the plot of the possible parameters between the right interval

"""# Surface of possible parameter space to explore
p_range = np.linspace(0, 1, 100)
mu_range = np.linspace(-np.log(8), np.log(8), 100)
sigma_range = np.linspace(0, np.log(16), 100)  # Adjusted range for sigma_1
p_grid, mu_grid, sigma_grid = np.meshgrid(p_range, mu_range, sigma_range)

# Plot the volume
ax.plot_trisurf(p_grid.ravel(), mu_grid.ravel(), sigma_grid.ravel(), linewidth=0, color='gray', alpha=0.2)"""

# Format and legend
ax.set_xlabel('$p$')
ax.set_ylabel('$\mu_1$')
ax.set_zlabel('$\sigma_1$')
if gibbs == 0:
    ax.set_title('Parameter space $(p, \mu_{1}, \sigma_{1})$ for MH algorithm')
else:
    ax.set_title('Parameter space $(p, \mu_{1}, \sigma_{1})$ for MH algorithm with Gibbs sampler')
cbar = plt.colorbar(sc)
cbar.set_label('Iterations')
plt.show()


## 2D animation graphs for sigma_1 in function of mu_1 at fixed p
"""
list_fixed_p = np.linspace(0.1, 0.9, 10)
list_fixed_mu_1 = np.linspace(-np.log(8), np.log(8), 10)
list_fixed_sigma_1 = np.linspace(0, np.log(16), 10)

# Set fixed axis limits
xmin = -np.log(8)
xmax = np.log(8)
ymin = 0
ymax = np.log(16)

# Create a list to store individual frames
frames = []

# Define the confidence level for the ellipses (e.g., 95%)
confidence_level = 0.95

# Create a colormap ranging from light red to dark red
cmap = plt.cm.get_cmap('Reds')

# Create the figure and axes outside the update_plot() function
fig, ax = plt.subplots()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

def update_plot(i):
    # Clear the previous plot
    ax.clear()

    fixed_p = list_fixed_p[i]
    fixed_mu_1 = 'None'
    fixed_sigma_1 = 'None'
    sample_MH = metropolis_hastings(normalized_treatment, x_values, num_samples, initial_parameters, var_m, gibbs, fixed_p, fixed_mu_1, fixed_sigma_1)
    liste_mu_1 = sample_MH[:, 1]
    liste_sigma_1 = sample_MH[:, 2]

    # Create the scatter plot for the new values of liste_mu_1 and liste_sigma_1
    scatter = ax.scatter(liste_mu_1, liste_sigma_1, c=np.linspace(0, 1, len(liste_mu_1)), marker='o', alpha=0.7)

    # Set the title for each graph with the value of fixed_p
    ax.set_title(f'Fixed $p$ = {fixed_p:.2f}')
    ax.set_xlabel('$\mu_{1}$')
    ax.set_ylabel('$\sigma_{1}$')

    # Calculate and display the covariance matrix and mean
    cov_matrix = np.cov(liste_mu_1, liste_sigma_1)
    mean_mu_1 = np.mean(liste_mu_1)
    mean_sigma_1 = np.mean(liste_sigma_1)

    # Draw ellipses representing the covariance matrix with the specified confidence level
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))

    # Calculate multiple ellipses for the covariance matrix
    num_ellipses = 4
    for k in range(1, num_ellipses + 1):
        ellipse_color = cmap(k / num_ellipses)  # Get the color from the colormap
        ellipse = Ellipse(
            xy=(mean_mu_1, mean_sigma_1),
            width=2 * np.sqrt(eigenvalues[0]) * k,
            height=2 * np.sqrt(eigenvalues[1]) * k,
            angle=angle,
            edgecolor=ellipse_color,
            facecolor='none'
        )
        ax.add_patch(ellipse)

    # Save the current frame as an image
    frames.append([ax])


# Create the animation
ani = FuncAnimation(fig, update_plot, frames=len(list_fixed_p)-2, interval=1000)

# Display the animation in the Jupyter Notebook
plt.show()

# Save the animation as a GIF using Pillow
ani.save('parameter_space_animation_mu_1_sigma_1.gif', writer='pillow')

"""

## 2D animation graphs for p in function of sigma_1 at fixed mu_1
"""
list_fixed_p = np.linspace(0.1, 0.9, 10)
list_fixed_mu_1 = np.linspace(-np.log(8), np.log(8), 10)
list_fixed_sigma_1 = np.linspace(0, np.log(16), 10)

# Set fixed axis limits
xmin = 0
xmax = np.log(16)
ymin = 0
ymax = 1

# Create a list to store individual frames
frames = []

# Define the confidence level for the ellipses (e.g., 95%)
confidence_level = 0.95

# Create a colormap ranging from light red to dark red
cmap = plt.cm.get_cmap('Reds')

# Create the figure and axes outside the update_plot() function
fig, ax = plt.subplots()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

def update_plot(i):
    # Clear the previous plot
    ax.clear()

    fixed_p = 'None'
    fixed_mu_1 = list_fixed_mu_1[i]
    fixed_sigma_1 = 'None'
    sample_MH = metropolis_hastings(normalized_treatment, x_values, num_samples, initial_parameters, var_m, gibbs, fixed_p, fixed_mu_1, fixed_sigma_1)
    liste_p = sample_MH[:, 0]
    liste_sigma_1 = sample_MH[:, 2]

    # Create the scatter plot for the new values of liste_mu_1 and liste_sigma_1
    scatter = ax.scatter(liste_sigma_1, liste_p, c=np.linspace(0, 1, len(liste_p)), marker='o', alpha=0.7)

    # Set the title for each graph with the value of fixed_p
    ax.set_title(f'Fixed $\\mu_{1}$ = {fixed_mu_1:.2f}')
    ax.set_xlabel('$\sigma_{1}$')
    ax.set_ylabel('$p$')

    # Calculate the covariance matrix and mean
    cov_matrix = np.cov(liste_sigma_1, liste_p)
    mean_p = np.mean(liste_p)
    mean_sigma_1 = np.mean(liste_sigma_1)

    # Draw ellipses representing the covariance matrix with the specified confidence level
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))

    # Calculate multiple ellipses for the covariance matrix
    num_ellipses = 4
    for k in range(1, num_ellipses + 1):
        ellipse_color = cmap(k / num_ellipses)  # Get the color from the colormap
        ellipse = Ellipse(
            xy=(mean_sigma_1, mean_p),
            width=2 * np.sqrt(eigenvalues[0]) * k,
            height=2 * np.sqrt(eigenvalues[1]) * k,
            angle=angle,
            edgecolor=ellipse_color,
            facecolor='none'
        )
        ax.add_patch(ellipse)

    # Save the current frame as an image
    frames.append([ax])


# Create the animation
ani = FuncAnimation(fig, update_plot, frames=len(list_fixed_mu_1)-2, interval=1000)

# Display the animation in the Jupyter Notebook
plt.show()

# Save the animation as a GIF using Pillow
ani.save('parameter_space_animation_p_sigma_1.gif', writer='pillow')
"""

## Animation of parameter space in 3D graph
"""
parametres = samples
iterations = np.arange(num_samples)

# Séparation des paramètres en trois listes pour chaque dimension
liste_p = parametres[:, 0]
liste_mu_1 = parametres[:, 1]
liste_sigma_1 = parametres[:, 2]

# Création de la figure 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Représentation des points en 3D
# sc = ax.scatter(liste_p, liste_mu_1, liste_sigma_1, c='black', marker='o', s=5, alpha=0.7)

# Dégradé de couleur selon les itérations (utilisation d'une carte de couleur "plasma" pour illustrer)
sc = ax.scatter(liste_p, liste_mu_1, liste_sigma_1, c=iterations, cmap='plasma', marker='o', s=10, alpha=0.7)

# Format
ax.set_xlabel('$p$')
ax.set_ylabel('$\mu_1$')
ax.set_zlabel('$\sigma_1$')
ax.set_title('Parameter space')
cbar = plt.colorbar(sc)
cbar.set_label('Iterations')

# Function to update the scatter plot for each iteration
def update_plot(i):
    # Clear the points for the previous iteration
    sc._offsets3d = (liste_p[:i], liste_mu_1[:i], liste_sigma_1[:i])
    # Update the colors for the previous iteration points
    sc.set_array(iterations[:i])
    # Set the title for the current iteration
    ax.set_title(f'Parameter space (Iteration {i})')

# Create the animation
ani = FuncAnimation(fig, update_plot, frames=num_samples, interval=100)
plt.show()

# Save the animation as a GIF using Pillow
ani.save('parameter_space_animation.gif', writer='pillow')"""
