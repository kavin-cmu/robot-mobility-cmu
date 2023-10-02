close all;
clear all;
clc;


% Model Parameters
lf = 1.50;
lr = 1.50;
A = [lf, lr]';

%% Part A

% Control Inputs
v = 2.50;
df = linspace(-pi/4,pi/4,10);

% Simultaion Time
dt = 0.05; Tf = 5.0;
ts = linspace(0, Tf, Tf/dt);

% Logging
X_log = zeros(3,length(ts),length(df));

X = [0,0,0]';

for i = 1:length(df)
    X = [0,0,0]';
    for j = 1:length(ts) 
        X = simulate_step(@pepyKBM,A,X,[v,df(i)],dt);
        X_log(:,j,i) = X';
    end
end

figure;
title('Path for Various Steering Inputs')
hold on;
grid on;
axis equal;

for i = 1:length(df)
    plot(X_log(1,:,i),X_log(2,:,i),'DisplayName',strcat('delta=',num2str((df(i)*180/pi))))
end

xlabel("X [m]")
ylabel("Y [m]")
title('Path for Various Steering Inputs')
legend('show')

figure;
hold on;
grid on;
% axis equal;

for i = 1:length(df)
    plot(ts, X_log(1,:,i),'DisplayName',strcat('delta=',num2str((df(i)*180/pi))))
end

xlabel("t [s]")
ylabel("X [m]")
legend('show')
title('X Trajectory')

figure;
hold on;
grid on;
% axis equal;

for i = 1:length(df)
    plot(ts, X_log(2,:,i),'DisplayName',strcat('delta=',num2str((df(i)*180/pi))))
end

xlabel("t [s]")
ylabel("Y [m]")
legend('show')
title('Y Trajectory')

figure;
hold on;
grid on;

% axis equal;

for i = 1:length(df)
    plot(ts, 180/pi*X_log(3,:,i),'DisplayName',strcat('delta=',num2str((df(i)*180/pi))))
end

xlabel("t [s]")
ylabel("\Phi [deg]")
legend('show')
title('\Phi Trajectory')

%% Part B


% Simulation Time
dt = 0.05; Tf = 5.0;
ts = linspace(0, Tf, Tf/dt);

% Control Inputs
v = 2.50;
df_amp = [pi/8, pi/6, pi/4];
df = zeros(3,length(ts));
df_freq = 0.5;

for i=1:length(df_amp)
    df(i,:) = df_amp(i)*sin(2*pi*df_freq*ts);
end

% Logging
X_log = zeros(3,length(ts), length(df));

X = [0,0,0]';

for i = 1:length(df_amp)
    X = [0,0,0]';
    for j = 1:length(ts) 
        X = simulate_step(@pepyKBM,A,X,[v,df(i,j)],dt);
        X_log(:,j,i) = X';
    end
end

figure;
hold on;
grid on;
% axis equal;

for i = 1:length(df_amp)
    plot(X_log(1,:,i),X_log(2,:,i),'DisplayName',strcat('delta_{amp} =',num2str((df_amp(i)*180/pi))))
end

xlabel("X [m]")
ylabel("Y [m]")
legend('show')
title('Path for Sinusoidal Steering Input')

figure;
hold on;
grid on;
% axis equal;

for i = 1:length(df_amp)
    plot(ts,X_log(3,:,i),'DisplayName',strcat('delta_amp=',num2str((df_amp(i)*180/pi))))
end

ylabel("\delta_f [rad]")
xlabel("t [s]")
legend('show')
title('Steering Input (\delta_f) ')

% figure;
% plot(ts,X)

%% Part C


% Simulation Time
dt = 0.05; Tf = 5.0;
ts = linspace(0, Tf, Tf/dt);

% Control Inputs
v = 2.50;
df_amp = [pi/8, pi/6, pi/4];
df = zeros(3,length(ts));
df_freq = 0.5;

for i=1:length(df_amp)
    df(i,:) = df_amp(i)*square(2*pi*df_freq*ts);
end

% Logging variable
X_log = zeros(3,length(ts), length(df));

X = [0,0,0]';

for i = 1:length(df_amp)
    X = [0,0,0]';
    for j = 1:length(ts) 
        X = simulate_step(@pepyKBM,A,X,[v,df(i,j)],dt);
        X_log(:,j,i) = X';
    end
end

figure;
hold on;
grid on;
axis equal;

for i = 1:length(df_amp)
    plot(X_log(1,:,i),X_log(2,:,i),'DisplayName',strcat('delta_{amp}=',num2str((df_amp(i)*180/pi))))
end

xlabel("X [m]")
ylabel("Y [m]")
legend('show')
title('Path for Square-waved Steering Input')

figure;
hold on;
grid on;
% axis equal;

for i = 1:length(df_amp)
    plot(ts,df(i,:),'DisplayName',strcat('delta_{amp}=',num2str((df_amp(i)*180/pi))))
end

ylabel("\delta_f [rad]")
xlabel("t [s]")
legend('show')
title('Steering Input (\delta_f) ')
%%
figure

for i = 1:length(df_amp)
    plot(ts,df(i,:))
    hold on;
    plot(ts,X_log(3,:,i),'DisplayName',strcat('delta_{amp}=',num2str((df_amp(i)*180/pi))))
end


