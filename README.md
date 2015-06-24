# entroloss 
Estimates entropy-loss in replica-averaged models using input from one (or more) canonical (unrestrained/unbiased) simulation(s) 

## usage
### input 
 * text-file with tabulated list of back-predicted data (dimensionality M) in each frame (J frames) from unbiased simulation(s) (will be handled as MxJ numpy array)
 * text-file with target mean (J)
 * residual arguments are integer values of the number of replicas to emulate N+1. e.g. the list  1 3 7 will evaluate the mean entropy-loss for 2 4 and 8 replicas.
### output
 * standard output with N+1 and entropy-loss 


For licensing see LICENSE
