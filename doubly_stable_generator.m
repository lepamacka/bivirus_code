%This script makes use of perturb_doublystable.m to generate examples 
%of bi-virus systems where both viruses have stable boundary equilibria. 
%The script requires a positive square matrix A, the spread matrix of 
%virus 1, from which B, the spread matrix of virus 2, is generated. 
%The script also approximates x and y, the corresponding boundary equilibria.

%Note that the returned matrices, boundary equilibria and eigenvalues 
%are subject to numerical errors. To increase confidence in a pair A and B, 
%choose a more negative eigmax_TOL. With more negative eigmax_TOL a larger 
%deltax_scale might be necessary to find a solution. Note that an excessively 
%negative eigmax_TOL can slow down or stall the program. 

clear all

eigmax_TOL = -0.00000001; %Limits the largest real part of any eigenvalue 
                          %of the Jacobians at the boundary equilibria.
                          %Lower limits can slow down or stall the program.

deltax_scale = 0.01;      %Determines the maximum scale of the difference 
                          %between x and y as a function of x.

n=5;                      %Number of nodes in the system. The program is 
                          %slower with larger n.
I = eye(n);

A = zeros(n,n);           %The spread matrix of virus 1.
while max(eig(A)) < 1
    A = 2*rand(n,n);      %Here A is randomized, it can also be defined freely.
end

%The following function call finds B, the spread matrix of virus 2, as well
%as approximations of the boundary equilibria x and y.
[x, y, B] = doubly_stable_func(A,deltax_scale,eigmax_TOL); 

A %Here the resulting pair of spread matrices are printed in output.
B

%The following line approximates the eigenvalues determining the stability 
%of the boundary equilibria x and y, and prints these eigenvalues in output.
stability_eigenvalues = [eigs(-I + (I-diag(y))*A, 1, 'largestreal') eigs(-I + (I-diag(x))*B, 1, 'largestreal')] 

