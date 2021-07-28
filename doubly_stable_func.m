%This function takes the first input, a spread matrix A, finds its boundary 
%equilibrium, x, and generates a perturbation of this equilibrium, deltax.
%It then uses deltax to find a spread matrix B with a boundary equilibrium y. 
%The end result is a bi-virus system defined by A and B with two stable 
%boundary equilibria, x and y.

%The second input is deltax_scale, defining the relative scale of the
%perturbation of x. Smaller deltax_scale than 1 will limit the potential
%magnitude of the stability-determining eigenvalues, but can also speed up
%the running time of this function.

%The third input is eigmax_TOL, which is an upper limit on the stability-
%determining eigenvalues. More negative eigmax_TOL increase the running
%time of this function and might cause it to stall.

function [x, y, B] = doubly_stable_func(A,deltax_scale,eigmax_TOL)

n=length(A);
I = eye(n);

x_zero = rand(n,1);
tspan = [0, 1e15];
[T,sol] = ode23s(@(t,p)(-I + (I-diag(p))*A)*p,tspan,x_zero);
x = sol(end, :)'; %This is an approximation of the boundary equilibrium of A.
[V1, D1, W1] = eig((I-diag(x))*A);
[maxeig1, ind1] = max(diag(D1));
u = W1(:,ind1)'/(W1(:,ind1)'*x); %This is the left Perron eigenvector of (I-diag(x))*A.

y = zeros(n,1);
B = zeros(n,n);
while eigs(-I + (I-diag(y))*A, 1, 'largestreal')>eigmax_TOL || eigs(-I + (I-diag(x))*B, 1, 'largestreal')>eigmax_TOL
zeta = 1-2*rand(n-1,1);
zeta(n) = -zeta' * x(1:n-1); %The vector zeta is perpendicular to x.
epsilon = 1;
B_p = A + epsilon * ones(n,1)*zeta'; %This is a spread matrix with x as its boundary equilibrium.
while ~all(B_p>0,'all') %This while loop ensures that B_p is a positive matrix.
    epsilon = epsilon*0.8;
    B_p = A + epsilon * ones(n,1)*zeta';
end

[V2, D2, W2] = eig((I-diag(x))*B_p);
[maxeig2, ind2] = max(diag(D2));
v = W2(:,ind2)'/(W2(:,ind2)'*x); %This is the left Perron eigenvector of (I-diag(x))*B_p.

I_neg_X_inv = diag((1-x).^(-1));
deltax = zeros(n,1); %This vector will be used to perturb B_p and is 
                     %randomly generated until it fulfils the inequalities 
                     %in the following while loop. The inequalities are 
                     %based on u and v, and will ensure double stability of 
                     %the final bi-virus system as long as deltax is 
                     %sufficiently small.
while (u*I_neg_X_inv*diag(x)*deltax <= 0 || v*I_neg_X_inv*diag(x)*deltax >= 0)
    deltax = deltax_scale*min([x 1-x],[],2).*(1-2*rand(n,1));
end

deltaB = zeros(n,n); %This matrix is generated using deltax in the following 
                     %steps, and will be used to perturb B_p.
targetvec = ((I-diag(x))^(-2)-B_p)*deltax;
rands2 = [zeros(n,1) sort(rand(n,n-1),2) ones(n,1)];
sum_one_vecs2 = rands2(:,2:(n+1))-rands2(:,1:n);
for i=1:n
    for j=1:n
        deltaB(i,j) = sum_one_vecs2(i,j)*x(i)/(targetvec(j));
    end
end

B = B_p+deltaB; %B is the spread matrix of virus 2.
while ~all(B>0,'all') %This while loop ensures that B is a positive matrix.
    deltaB = deltaB*0.8;
    B = B_p+deltaB;
end

[T,sol] = ode23s(@(t,p)(-I + (I-diag(p))*B)*p,tspan,x_zero);
y = sol(end, :)'; %This is an approximation of the boundary equilibrium of B.

ind = 0; %The following while loop checks that x and y are stable boundary 
         %equilibria in the system defined by A and B, and if not, deltaB
         %is decreased until the condition is met. The index ind is used to
         %break the while loop if it does not complete, which restarts the
         %search for a suitable B.
while ((eigs(-I + (I-diag(y))*A, 1, 'largestreal')>eigmax_TOL || eigs(-I + (I-diag(x))*B, 1, 'largestreal')>eigmax_TOL) && (~all(B>0,'all')) && (ind<1e3))
    deltaB = deltaB*0.8;
    B = B_p+deltaB;
    [T,sol] = ode23s(@(t,p)(-I + (I-diag(p))*B)*p,tspan,x_zero);
    y = sol(end, :)';
end
end
end

