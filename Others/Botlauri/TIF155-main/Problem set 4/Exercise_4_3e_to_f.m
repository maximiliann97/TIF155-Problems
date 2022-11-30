%% e) Calculation of the Lyapunov exponents.
% Parameters.
clearvars
tic
a = 1.4; b = 0.3;
nMax = 10^6;
ns = linspace(1, nMax, nMax);

% Initialization.
x0 = 0.1; y0 = 0.1;
x = zeros(1, nMax); y = zeros(1, nMax);
x(1) = x0; y(1) = y0;

% Calculation of the trajectories. 
for n = 1:(nMax - 1)
    x(n + 1) = y(n) + 1 - a*x(n).^2;
    y(n + 1) = b*x(n);
end
Q = eye(2);
lambda = zeros(1, 2);
lambdaList = zeros(nMax, 2);

% Calculation of the eigenvalues.
for i = 1:nMax
    traj = [x(i), y(i)];
    J = [-2*a*traj(1), 1; b, 0];
    M = J*eye(2);
    [Q,R] = qr(M*Q);
    lambda(1) = lambda(1) + 1/nMax*log(abs(R(1,1)));
    lambda(2) = lambda(2) + 1/nMax*log(abs(R(2,2)));
    lambdaList(i, :) = nMax/i*lambda;
end

% f) Calculation of the Lyapunov exponents using the Kaplan-Yorke conjecture. 
DL = 1 - lambda(1)/lambda(2);

toc

%% Plot of the development of the Lyapunov exponents.
plot(log(ns), lambdaList)
title('Approximation of the Lyapunov exponents \lambda_1 and \lambda_2.')
