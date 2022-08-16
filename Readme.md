# Usage ReadMe

### Usage

Download the snake_repo_parallel folder

Add snaker_repo_parallel and all subfolders to the MATLAB Path

run run_synthesis.m

You will be prompted to enter 7, 15, or 23, corresponding to the number of to the number of the transition system you wish to compute.

If you wish to test and System Labeling or Find Controlled Invariant Set functions, simply uncomment them. WARNING: these functions have not been fully modified to support the parallelized System Loading function. 
 
#### Notes

1. Note that there is a small bug causing the parallelized transition system to be slightly different from the original iteratively built transition system

2. Only the System Loading function has been parallelized. More work is needed to parallelize the Find Controlled Invariant Set function.

3. IMPORTANT: run_sim.m, while intended to be a visual representation of the transition system and its corresponding CIS, has not been updated and is not compatible with the parallelized code.
