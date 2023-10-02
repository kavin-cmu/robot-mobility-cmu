
function Xn = simulate_step(dynamics, A, X, U,dt)
    %simulate_step - updates the system state for a single simulation step
    %
    % Syntax:  Xn = simulate_step(dynamics, A, X, U,dt)
    %
    % Inputs:
    %    dynamics - function handle to the state derivative function
    %    A - KBM model parameters
    %    X - current vehicle state
    %    U - control input vector
    %    dt - simulation timestep
    %
    % Outputs:
    %    Xn - updated state
    
     Xn = X + dynamics(A, X, U)*dt;
 end