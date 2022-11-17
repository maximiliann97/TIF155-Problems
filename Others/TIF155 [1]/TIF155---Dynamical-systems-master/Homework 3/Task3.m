clear
close all
sigma = 10;
b = 8/3;
r = 0.32;

threshold = 30000;
deltaT = 0.0001;
transientBuffer = 10000;
maxT = transientBuffer*deltaT + 1;
timeInterval = transientBuffer*deltaT:deltaT:maxT;
totIter = floor(transientBuffer*deltaT + threshold/deltaT);

x0 = 0.1;
y0 = 0.1;
z0 = 0.1;

% Step 1: Solve x' = f(x)
f = @(t, var) [sigma*(var(2)-var(1)); r*var(1)-var(2)-var(1)*var(3);var(1)*var(2)-b*var(3)];
[t,var] = ode45(f, timeInterval, [x0,y0,z0]);
maxIter = size(var,1);

% Step 2: Define Q as identity matrix, and lambdas as zero
Q = eye(3);
lambda = zeros(3, totIter);
lambdaSum = zeros(1, totIter);

timeStep = transientBuffer*deltaT;
iter = transientBuffer+1;
tempIter = 1;
while timeStep <= threshold
    if mod(iter, 10000) == 0
        disp(iter);
        lambdaFinSum = sum(lambda,2)/((iter-1-transientBuffer)*deltaT);
        disp(lambdaFinSum);
        disp(sum(lambdaFinSum));
        disp(timeStep);
    end
    
    if tempIter >= maxIter
        maxT = maxT + 1;
        timeInterval = maxT-1+deltaT:deltaT:maxT;
        
        % Step 1: Solve x' = f(x)
       % f = @(t, var) [sigma*(var(2)-var(1)); r*var(1)-var(2)-var(1)*var(3);var(1)*var(2)-b*var(3)];
        [t,var] = ode45(f, timeInterval, [var(tempIter-1,1),var(tempIter-1,2),var(tempIter-1,3)]);
        maxIter = size(var,1);
        tempIter = 1;
    end
    % Step 3: alculate new M
    J = [-sigma, sigma, 0; r-var(tempIter,3), -1, -var(tempIter,1); var(tempIter,2), var(tempIter,1), -b];
    M = eye(3)+J.*deltaT;
    
    % Step 4: 
    [Q,R] = qr(M*Q);
    
    % Step 5: Add diagonal components of R to lambda
    diagElements = diag(R);
    
    lambda(1, iter) = log(abs(R(1,1)));
    lambda(2, iter) = log(abs(R(2,2)));
    lambda(3, iter) = log(abs(R(3,3)));
    lambdaSum(iter) = sum(log(abs(diagElements)));
    
    timeStep = timeStep + deltaT;
    iter = iter+1;
    tempIter = tempIter + 1;
end
%lambdaFin1 = sum(lambda(1,:))/(timeStep-transientBuffer*deltaT);
%lambdaFin2 = sum(lambda(2,:))/(timeStep-transientBuffer*deltaT);
%lambdaFin3 = sum(lambda(3,:))/(timeStep-transientBuffer*deltaT);
%lambdaFinSum = sum(lambdaSum)/(timeStep-transientBuffer);
lambdaFinSum = sum(lambda,2)/((iter-1-transientBuffer)*deltaT);
lambdaArray(1,:) = cumsum(lambda(1,(transientBuffer+1):iter-1))./(((transientBuffer+1):iter-1)*deltaT);
lambdaArray(2,:) = cumsum(lambda(2,(transientBuffer+1):iter-1))./(((transientBuffer+1):iter-1)*deltaT);
lambdaArray(3,:) = cumsum(lambda(3,(transientBuffer+1):iter-1))./(((transientBuffer+1):iter-1)*deltaT);

plot3(var(:,1),var(:,2),var(:,3))
semilogx(transientBuffer+1:iter-1,lambdaArray)
sum(lambdaFinSum);