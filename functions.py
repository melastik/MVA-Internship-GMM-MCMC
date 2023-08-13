import numpy as np
from scipy.stats import norm
from tqdm import tqdm
import scipy
import scipy.stats
import csv
import subprocess
import rpy2.robjects as robjects

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import ArtistAnimation
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
from PIL import Image
from IPython.display import HTML

############################# Useful functions #############################

# Function to get .csv row into a list with values and p-value
def get_csv_line(file_path, line_number):
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file, delimiter=';')
        for i, row in enumerate(csv_reader):
            if i == line_number:
                converted_row = [row[0], row[1], int(row[2]), int(row[3]), int(row[4]), int(row[5]), int(row[6]), int(row[7]), int(row[8]), int(row[9]), int(row[10]), float(row[11])]
    return converted_row

# Function to get the row content of a csv file in a table
def data(csv_file,n):
    line_table = []
    for line_number in tqdm(range(1,n)): # 11622
        # print('Row number in .csv file: ', line_number) without tqdm
        c_row = get_csv_line(csv_file, line_number)
        line_table.append([c_row[2], c_row[3], c_row[4], c_row[5], c_row[6], c_row[7], c_row[8], c_row[9], c_row[10]])
    return line_table
    
# Function to sum mutants on a line i for a gene
def sum_mutants_gene(table_list,i):
    sum_mutants = np.sum(table_list[i])
    return sum_mutants

# Function to normalize cumulative line of a choosen treatment
def normalization(treatment):
    return treatment / np.sum(treatment)

# Function to compute the probability density of the Gaussian distribution
def gaussian_distribution(x, mu, sigma):
    return (1 / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

# Function to compute the proportion of values in the target histogram within the interval [0.5, 2]
def proportion_in_interval(target_hist, center_values):
    interval_mask = (center_values >= 0.5) & (center_values <= 2)
    return np.sum(target_hist[interval_mask])

# Function for the proposal distribution (gaussian mixture)
def gaussian_mixture(x, p, mu_0, sigma_0, mu_1, sigma_1):
    return p * (1 / np.sqrt(2 * np.pi * sigma_0 ** 2)) * np.exp(-((x - mu_0) ** 2) / (2 * sigma_0 ** 2)) + \
           (1 - p) * (1 / np.sqrt(2 * np.pi * sigma_1 ** 2)) * np.exp(-((x - mu_1) ** 2) / (2 * sigma_1 ** 2))

# Function for CDF of gaussian distribution of parameters mu and sigma
def norm_cdf(x,mu,sigma):
    # vec_cdf = norm.cdf(x_values)
    return norm.cdf((x-mu)/sigma)
    
# Function to compute proposed histogram []*9 (vecteur y proposé)
def y_prop(vec, p, mu_0, sigma_0, mu_1, sigma_1):
    vec_cdf_0 = norm_cdf(vec,mu_0,sigma_0) # gaussian 1
    vec_cdf_1 = norm_cdf(vec,mu_1,sigma_1) # gaussian 2
    end = len(vec_cdf_0)
    # print("vec_cdf_0", vec_cdf_0[1:end]-vec_cdf_0[0:end-1])
    return p*(vec_cdf_0[1:end]-vec_cdf_0[0:end-1])+(1-p)*(vec_cdf_1[1:end]-vec_cdf_1[0:end-1])

# MH with Gibbs sampler fonctionnel, fixed parameter and 3D plot
def metropolis_hastings(target_hist, x_values, num_samples, mu_0, sigma_0, initial_params, var_m, gibbs, fixed_p, fixed_mu_1, fixed_sigma_1): # fixed_p, fixed_mu_1, fixed_sigma_1
    
    # Initialization
    p_current, mu_1_current, sigma_1_current = initial_params
    accepted_samples = np.zeros((num_samples,3))
    y_current = y_prop(x_values, p_current, mu_0, sigma_0, mu_1_current, sigma_1_current)
        
    # With Gibbs sampler
    if gibbs == 1:
        for i in range(num_samples):
            p_proposed, mu_1_proposed,sigma_1_proposed = p_current,mu_1_current,sigma_1_current
            if i%3 == 0 :
                p_proposed = np.random.uniform(0, 1, 1)
            if i%3 == 1 :
                mu_1_proposed = np.random.uniform(-np.log(8), np.log(8), 1)
            if i%3==2 :
                sigma_1_proposed = np.random.uniform(0, np.log(16), 1)
            y_proposed = y_prop(x_values, p_proposed, mu_0, sigma_0, mu_1_proposed, sigma_1_proposed)

            # Compute the acceptance ratio
            # acceptance_ratio = (np.sum((target_hist - y_current)**2) - np.sum((target_hist - y_proposed)**2))/(sigma_m)**2
            acceptance_ratio = (np.sum((target_hist - y_current)**2) - np.sum((target_hist - y_proposed)**2))/var_m

            # Accept or reject the proposed parameters
            if np.log(np.random.rand()) < acceptance_ratio: # np.log(np.random.rand()) est nécessairement négatif
                p_current, mu_1_current, sigma_1_current = p_proposed, mu_1_proposed, sigma_1_proposed

            y_current = y_prop(x_values, p_current, mu_0, sigma_0, mu_1_current, sigma_1_current)
            accepted_samples[i] = p_current, mu_1_current, sigma_1_current
            
        return np.array(accepted_samples)
        
    # Without Gibbs sampler
    if gibbs == 0:
        
        for i in range(num_samples):
            
            if fixed_p == 'None': # for non fixed parameters p
                p_proposed = np.random.uniform(0, 1, 1)
            else:
                p_proposed = fixed_p
                
            if fixed_mu_1 == 'None': # for non fixed parameters mu_1
                mu_1_proposed = np.random.uniform(-np.log(8), np.log(8), 1)
            else:
                mu_1_proposed = fixed_mu_1
                
            if fixed_sigma_1 == 'None': # for non fixed parameters sigma_1
                sigma_1_proposed = np.random.uniform(0, np.log(16), 1)
            else:
                mu_1_proposed = fixed_sigma_1
            
            y_proposed = y_prop(x_values, p_proposed, mu_0, sigma_0, mu_1_proposed, sigma_1_proposed)

            # Compute the acceptance ratio
            # acceptance_ratio = (np.sum((target_hist - y_current)**2) - np.sum((target_hist - y_proposed)**2))/(sigma_m)**2
            acceptance_ratio = (np.sum((target_hist - y_current)**2) - np.sum((target_hist - y_proposed)**2))/var_m

            # Accept or reject the proposed parameters
            if np.log(np.random.rand()) < acceptance_ratio: # np.log(np.random.rand()) < 0
                p_current, mu_1_current, sigma_1_current = p_proposed, mu_1_proposed, sigma_1_proposed

            y_current = y_prop(x_values, p_current, mu_0, sigma_0, mu_1_current, sigma_1_current)
            
            accepted_samples[i] = p_current, mu_1_current, sigma_1_current
        
    return np.array(accepted_samples)

# Function to simulate according to the estimated gaussian mixture model and from the data table_list
def simulation_mutants(table_list, n_row, mu_0, sigma_0, estimated_mu_1, estimated_sigma_1, estimated_p):
    sum_mutants = [0]*n_row
    z_list = []
    simulated_table_list = []
    for line in range(n_row-1):
        sum_mutants[line] = sum_mutants_gene(table_list, line)
        
        Z = np.random.binomial(1, 1-estimated_p)
        
        if Z == 0: # If Z is 0, the sample comes from the first Gaussian (Gaussian 1)
            z_list.append(Z)
            sample = np.random.normal(mu_0, sigma_0, sum_mutants[line])
            simulated_table_list.append(sorted(sample.tolist()))
        else: # If Z is 1, the sample comes from the second Gaussian (Gaussian 2)
            z_list.append(Z)
            sample = np.random.normal(estimated_mu_1, estimated_sigma_1, sum_mutants[line])
            simulated_table_list.append(sorted(sample.tolist()))
            
    return z_list, simulated_table_list

# Function to sample the mutants form the simulation in intervals
def sampling_mutants(x_values, simulated_table_list):
    # Initialize a list of lists to store the counts
    counts_list = [[0] * (len(x_values) - 1) for _ in range(len(simulated_table_list))]

    # Iterate through each sublist in simulated_table_list
    for i, sublist in enumerate(simulated_table_list):
        # Iterate through each value in the sublist
        for value in sublist:
            # Find the interval to which the value belongs
            interval_index = None
            for j in range(len(x_values) - 1):
                if x_values[j] <= value < x_values[j + 1]:
                    interval_index = j
                    break

            # If the value does not belong to any interval, skip it
            if interval_index is not None:
                # Increment the count of the corresponding interval for this sublist
                counts_list[i][interval_index] += 1

    return counts_list

# Function to compute the cumulative simulated histogram
def simulated_cml_hist(n_row, simulated_treatment_data):
    # Initialisation de la liste cml_hist avec des zéros
    cml_hist = [0] * len(simulated_treatment_data[0])
    
    for k in range(n_row-1):
        for i in range(len(simulated_treatment_data[0])):
            # Ajouter la valeur actuelle à la somme correspondante dans cml_hist
            cml_hist[i] += simulated_treatment_data[k][i]
            
    return cml_hist

        
