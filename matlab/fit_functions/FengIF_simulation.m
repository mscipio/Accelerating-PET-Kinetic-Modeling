%|=================================================================================
%|   THEORETICAL MODEL OF INPUT FUNCTION USED FOR SIMULATION
%|
%|   INPUTS:
%|       time:           time vector [time_points x 1]
%|       inputFunParam:  input function model parameters to simulate
%|                      (delay, A1, A2, A3, lambda1, lambda2, lambda3) [1 x 7]
%|
%|   OUTPUTS:
%|       Cp:       simulated AIF time curve [time_points x 1]
%|       ppCp:     spline approximation of simulated AIF time curve [time_points x 1]
%|       ifParams: input function model parameters to simulate
%|                 (delay, A1, A2, A3, lambda1, lambda2, lambda3) [1 x 7]
%|                 this takes care of missing values in the input vector
%|  Last revision:
%|  22 May 2018
%|  Michele Scipioni, Univeristy of Pisa
%|
%|=================================================================================



function [Cp, ppCp,ifParams] = FengIF_simulation(time, inputFunParam)

% With time(r,1) containing the time points in minutes
% ifParams: [delay, A1, A2, A3, -lambda1, -lambda2, -lambda3]
% Use togheter with FengInput.m

[Cp, ifParams] = FengInput([2  inputFunParam],time/60);
ppCp = spline(time,FengInput(2,time));

