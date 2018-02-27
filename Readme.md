# Readme
This folder contains the files and codes for 2 state markov switching Current Account Deficit analysis. 

To run the code copy the folder to your local and change the first line in main.R to set working directory. 
The code is written in R 3.4 64bit and requires the following packages:

## Requirements
- ggplot2
- gridExtra

## Files
- main.R: This is the main file. The results can be replicated by running main.R
- twostate_DSM,twostate_DM,twostate_DS,twostate_D: These files apply the 2 state MS-ARFIMA model to the given data, CAD_quarterly_BOP.rda. The models executed in the files allow (d,sigma,mu), (d,mu), (d,sigma), (d) to switch state respectively.
  Inputs:
    1. filename: Name of the data to be analyzed (an .rda file).
  Outputs:
    1. dataname2_stateD**_ress.csv: Parameter estimates, likelihoods and optimization reports for the fitted model. The results are saved to output folder.
- utils.R: Contains various functions to help parse the results:
  1. totable: Converts list of results to table
  2. correctRes: # Reformats the table of estimates in "dH,dL,PHH,PLL, ..." order, instead of "dL,dH,PLL,PHH, ..." order.
  3. report: reports statistics of the outputs. See the paper
  4. processdat: Returns a given series by excluding the NA's at the beginning and the end.
  5. deprocessds: According to what series they belong to, adds NA's at the beginning and the end of the estimated d series.

- plots.R: Plots the paths of estimated d series for 4 type of models (DSM,DM,DS,D). The plots are saved to pairplots folder.
