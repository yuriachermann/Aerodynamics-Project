clear all
close all
clc

%% Input data

% Flight conditions
flycond.rho = 1.225;
flycond.Uinf = 60;

% Wing characteristics
wing.AR = 12;
wing.S = 16.4;
wing.b = sqrt(wing.AR * wing.S);
wing.lambda = 0.56;
wing.rchord = 2 * wing.S / (wing.b * (1 + wing.lambda));
wing.tchord = wing.rchord * wing.lambda;
wing.stations = [-wing.b/2 0 wing.b/2];
wing.chords = [wing.tchord wing.rchord wing.tchord];

wing.airfoil.alpha_L_0 = -1.213;
wing.airfoil.clalpha = 6.13;
wing.alpha = 5;

%% Solve

N = 30;
A = LLT_solve(wing, flycond, N);

%% Post-processing

n = 100;
LLT_plot_results(A, flycond, wing, n);