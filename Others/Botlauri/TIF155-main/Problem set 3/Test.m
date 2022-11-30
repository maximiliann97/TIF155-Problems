% Test.
tic
sigma = 16; b = 5; r = 350;
Tmax = 1000; Tmin = 1; dt = 0.0001;
ts = linspace(Tmin, Tmax, (Tmax - Tmin)/dt);
[t,y] = ode45(@ode, ts, [0.1; 0.1; 0.1]);
Q = eye(3);
lambda = zeros(1, 3);
lambdaList = zeros(length(ts), 3);
for i = 1:length(ts)
    traj = y(i, :);
    J = [-sigma, sigma, 0; r-traj(3), -1, -traj(1); traj(2), traj(1), -b];
    M = eye(3) + J*dt;
    [Q,R] = qr(M*Q);
    lambda(1) = lambda(1) + 1/Tmax*log(abs(R(1,1)));
    lambda(2) = lambda(2) + 1/Tmax*log(abs(R(2,2)));
    lambda(3) = lambda(3) + 1/Tmax*log(abs(R(3,3)));
    lambdaList(i, :) = Tmax/(i*dt)*lambda;
end
plot(log(t), lambdaList)
toc

% The ODE.
function dydt = ode(t, y)
    sigma = 10; b = 8/3; r = 28;
    dydt = [sigma*(y(2) - y(1)); r*y(1)-y(2)-y(1)*y(3); y(1)*y(2)-b*y(3)];
end
