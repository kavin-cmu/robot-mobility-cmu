function Xdot = pepyKBM(A,X,U)
    %pepyKBM - calculates the state derivative given the model, current state and the control inputs
    %
    % Syntax:  Xdot = pepyKBM(A,X,U)
    %
    % Inputs:
    %    A - KBM model parameters
    %    X - current vehicle state
    %    U - control input vector
    %
    % Outputs:
    %    Xdot - current state derivative

    Xdot = [0,0,0]';
    Xdot(1) =  U(1)*cos(X(3));
    Xdot(2) =  U(1)*sin(X(3));
    Xdot(3) =  U(1)*tan(U(2))/(A(1)+A(2));
end