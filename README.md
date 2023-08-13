# Identifying the genetic determinants of chemical toxicity in microalga Chlamydomonas reinhardtii thanks to statistical analysis - Gaussian mixture models estimation 

## Abstract 
Microalgae such as *Chlamydomonas reinhardtii* offer significant potential for producing valuable compounds in industry (lipids, pigments, ...). This research project investigates genetic factors behind chemical toxicity responses in Chlamydomonas reinhardtii, focusing on herbicides atrazin, diuron and paraquat. Through experimental tests conducted on mutants of \textit{Chlamydomonas reinhardtii}, we studied the microalga's genome to discern genes exhibiting sensitivity or tolerance to these specific herbicides, determined by evaluating the mutant phenotype ratio. This ratio is an indicator of treatment impact. \\
The study is based on Gaussian mixture models (GMM) for statistical analysis, enabling the identification of genes with significant treatment responses. The parameters of GMM were estimated thanks to Metropolis-Hasting algorithm (MH). Various simulations were performed then to show the impact of sampling mutant phenotype ratio in experimental process and how it impacts the statistical analysis. 

## Content of the repository
- `functions.py` is a Python script containing the required functions to run the project.
- `main.py` is a Python script used to set parameters of functions defined in functions.py and plot results. 
- `Fisher-exact-test.R` is a R script used to run Fisher's exact test on simulated data from main.py and generate a .csv file that contains the p-values for a given treatment. 
- `data_example.csv` is a table containing fictive data to run Python script
- `requirements.txt` is a text file that contains all the versions of Python librairies used to run the whole project.

## Dependencies
**Hardware and working setup**
- macOS Monterey 12.6.1
- Anaconda 23.1.0 with Python 3.9.12

**Python librairies**
The librairy list and their versions are listed in `requirements.txt`. 
- Numpy 1.21.5
- Scipy 1.7.3
- tqdm 4.64.1
- Scipy 1.7.3
- csv 
- Subprocess
- Rpy2 3.5.13
- Matplotlib 3.5.3
- Matplotlib-inline 0.1.6
- Mpl_toolkits.mplot3d 
- Pillow 9.2.0
- iPython 8.6.0


## How to run the code ?
If you want to create an Anaconda environment to run the project, you need to install Anaconda with the following link https://www.anaconda.com/download/
Then, open Terminal and run the following lines to create a conda environment for this project :

```
conda create --name MVA-Internship-GMM-MCMC-env python=3.9.12
conda activate MVA-Internship-GMM-MCMC-env
```

Clone this repository to run the code: 
'git clone https://github.com/melastik/MVA-Internship-GMM-MCMC'

Go in the path you clone the repository:

`cd /Users/.../MVA-Internship-GMM-MCMC/`

Then you can install all the librairies required in the project using the following command line:

`pip install -r requirements.txt`

To run a demo:

`python main.py`


