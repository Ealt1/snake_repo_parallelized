%% Add Paths
%Citation info:
%https://www.mathworks.com/matlabcentral/fileexchange/60135-parfor-progress-monitor-progress-bar-v3

clear,clc


addpath('lib_abstr')
addpath('lib_cinv_synthesis')
addpath('lib_plot')
addpath('sys1')
addpath('ndSparse_G4_2021_03_16')

%% System Loading

%create a pool of workers
gcp

tic
%Create the system (dynamics)
while 1
    myinput = input("How many nodes would you like to process?\n Please enter either 7, 15, or 23, then hit Enter")
    if myinput == 7
        sys = node_7();
        break
    elseif myinput == 15
        sys = node_15();
        break
    elseif myinput == 23
        sys = node_23();
        break
    end
end
toc

%% System Labeling
% %Label system nodes as safe or unsafe
% tic
% X_safe = label_sys(sys); 
% toc

%% Find Controlled Invariant Set
%Find the controlled invariant set within the safe set
%Returns C (controlled invariant set), and K (the controller associated
%with set C)

% tic
% [C,K] = findControlledInvariantMulti_cell(X_safe, sys.tau);
% toc
