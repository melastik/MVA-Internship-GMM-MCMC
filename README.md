# Identifying the genetic determinants of chemical toxicity in microalga Chlamydomonas reinhardtii thanks to statistical analysis - Gaussian mixture models estimation 

## Content of the repository
- 'functions.py' is a Python script containing the required functions to run the project.
- 'main.py' is a Python script used to set parameters of functions defined in functions.py and plot results. 
- 'Fisher-exact-test.R' is a R script used to run Fisher's exact test on simulated data from main.py and generate a .csv file that contains the p-values for a given treatment. 
- 'data_example.csv' is a table containing fictive data to run Python script
- 'requirements.txt' is a text file that contains all the versions of Python librairies used to run the whole project.

## Requirements to run the code** 
The code requires Python 3.9.12. 

If you want to create an Anaconda environment to run the project, you need to install Anaconda with the following link https://www.anaconda.com/download/
Then, open Terminal and run the following lines to create a conda environment for this project :

'''
conda create --name MVA-Internship-GMM-MCMC-env python=3.9.12
conda activate MVA-Internship-GMM-MCMC-env
'''

Clone this repository to run the code: 
'git clone https://github.com/melastik/MVA-Internship-GMM-MCMC'

'cd /USers/MelaniePietri/MVA-Internship-GMM-MCMC/'

Then you can install all the librairies required in the project using the following command line:

'pip install -r requirements.txt'

To run a demo:

'python main.py'


