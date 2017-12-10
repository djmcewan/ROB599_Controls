global gridN
gridN = 20;
tic
% Minimize the simulation time
time_min = @(x) x(1);
% The initial parameter guess; 1 second, gridN positions, gridN velocities,
% gridN accelerations
x0 = [1; linspace(0,1,gridN)'; linspace(0,1,gridN)'; ones(gridN, 1) * 5];
% No linear inequality or equality constraints
A = [];
b = [];
Aeq = [];
Beq = [];
% Lower bound the simulation time at zero seconds, and bound the
% accelerations between -10 and 30
lb = [0;    ones(gridN * 2, 1) * -Inf;  ones(gridN, 1) * -10];
ub = [Inf;  ones(gridN * 2, 1) * Inf;   ones(gridN, 1) * 30];
% Options for fmincon
options = optimoptions(@fmincon, 'TolFun', 0.00000001, 'MaxIter', 10000, ...
                       'MaxFunEvals', 100000, 'Display', 'iter', ...
                       'DiffMinChange', 0.001, 'Algorithm', 'sqp');
% Solve for the best simulation time + control input
optimal = fmincon(time_min, x0, A, b, Aeq, Beq, lb, ub, ...
              @double_integrator_constraints, options);
% Discretize the times
sim_time = optimal(1);
delta_time = sim_time / gridN;
times = 0 : delta_time : sim_time - delta_time;
% Get the state + accelerations (control inputs) out of the vector
positions = optimal(2             : 1 + gridN);
vels      = optimal(2 + gridN     : 1 + gridN * 2);
accs      = optimal(2 + gridN * 2 : end);
% Make the plots
figure();
plot(times, accs);
title('Control Input (Acceleration) vs Time');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
figure();
plot(times, vels);
title('Velocity vs Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
figure();
plot(times, positions);
title('Position vs Time');
xlabel('Time (s)');
ylabel('Position (m)');
fprintf('Finished in %f seconds', toc);

function [ c, ceq ] = double_integrator_constraints( x )
    global gridN
    % No nonlinear inequality constraint needed
    c = [];
    % Calculate the timestep
    sim_time = x(1);
    delta_time = sim_time / gridN;
    % Get the states / inputs out of the vector
    positions = x(2             : 1 + gridN);
    vels      = x(2 + gridN     : 1 + gridN * 2);
    accs      = x(2 + gridN * 2 : end);
    
    % Constrain initial position and velocity to be zero
    ceq = [positions(1); vels(1)];
    for i = 1 : length(positions) - 1
        % The state at the beginning of the time interval
        x_i = [positions(i); vels(i)];
        % What the state should be at the start of the next time interval
        x_n = [positions(i+1); vels(i+1)];
        % The time derivative of the state at the beginning of the time
        % interval
        xdot_i = [vels(i); accs(i)];
        % The time derivative of the state at the end of the time interval
        xdot_n = [vels(i+1); accs(i+1)];
        
        % The end state of the time interval calculated using quadrature
        xend = x_i + delta_time * (xdot_i + xdot_n) / 2;
        % Constrain the end state of the current time interval to be
        % equal to the starting state of the next time interval
        ceq = [ceq ; x_n - xend];
    end
    % Constrain end position to 1 and end velocity to 0
    ceq = [ceq ; positions(end) - 1; vels(end)];
end
