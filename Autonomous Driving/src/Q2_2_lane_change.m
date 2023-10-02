close all;
clear all;
clc;

% Model Parameters
m = 1573;
Iz = 2873;
lf = 1.10;
lr = 1.58;
Cf = 8e4;
Cr = 8e4;
Vx = 30;

% System Matrix
A = 2*[ 0,  1/2,                    0,                  0;
        0, -(Cf+Cr)/(m*Vx),         (Cf+Cr)/m,          (-Cf*lf+Cr*lr)/(m*Vx);
        0,  0,                      0,                  1/2;
        0, -(Cf*lf-Cr*lr)/(Iz*Vx),  (Cf*lf-Cr*lr)/Iz,   -(Cf*lf^2+Cr*lr^2)/(Iz*Vx)];
%Control Matrix
B1 = [0;
     2*Cf/m;
     0;
     2*Cf*lf/Iz;
    ];

%Feed-Forward Matrix
B2 = [0;
      -2*(Cf*lf-Cr*lr)/(m*Vx) - Vx;
      0;
      -2*(Cf*lf^2+Cr*lr^2)/(Iz*Vx)
     ];

%% Reference Path
strt_seg_len1_x = 5.0;
strt_seg_len2_x = 35.0;
slant_seg_len_x = 90; 
total_path_len = strt_seg_len2_x+strt_seg_len1_x+slant_seg_len_x;
slant_seg_len_y = 5.0;
str_seg1_y = -5.0;
str_seg2_y = 0.0;
slant_ang = atan2(slant_seg_len_y,slant_seg_len_x);
slant_seg_len = norm([slant_seg_len_x,slant_seg_len_y]);

%% Time Parametrization
dt = 0.01;
tf = (slant_seg_len + strt_seg_len2_x + strt_seg_len1_x)/Vx;
t_strt1 = strt_seg_len1_x/Vx;
t_slant = t_strt1+slant_seg_len/Vx;
ts1 = linspace(0,t_strt1,floor(t_strt1/dt));
ts2 = linspace(t_strt1+dt,t_slant,floor((t_slant-t_strt1)/dt));
ts3 = linspace(t_slant+dt,tf,floor((tf-t_slant)/dt));
ts = [ts1, ts2, ts3];

dphi_wind = 10*dt;

x_ref = zeros(1,length(ts));
y_ref = zeros(1,length(ts));
dy_ref = zeros(1,length(ts));
phi_ref = zeros(1,length(ts));
dphi_ref = zeros(1,length(ts));
X_ref = zeros(4,length(ts));
x_act = zeros(1,length(ts));

x_act(1) = 0;

x_ref(1:length(ts1)) = linspace(0,strt_seg_len1_x,length(ts1));
x_ref(length(ts1)+1:length(ts1)+length(ts2)) = strt_seg_len1_x + linspace(Vx*dt*cos(slant_ang),slant_seg_len_x, length(ts2));
x_ref(length(ts1)+length(ts2)+1:length(ts)) = linspace(Vx*dt+strt_seg_len1_x+slant_seg_len_x,strt_seg_len1_x+slant_seg_len_x+strt_seg_len2_x,length(ts3));

y_ref(1:length(ts1)) = str_seg1_y;
y_ref(length(ts1)+1:length(ts1)+length(ts2)) = linspace(Vx*dt*sin(slant_ang)+str_seg1_y,str_seg2_y, length(ts2));
y_ref(length(ts1)+length(ts2)+1:length(ts)) = str_seg2_y;

dy_ref(1:length(ts1)) = 0.0;
dy_ref(length(ts1)+1:length(ts1)+length(ts2)) = Vx*sin(slant_ang);
dy_ref(length(ts1)+length(ts2)+1:length(ts)) = 0.0;

dphi1_start = floor((t_strt1-dphi_wind)/dt);
dphi1_end = floor((t_strt1)/dt);
dphi2_start = floor((t_slant-dphi_wind)/dt);
dphi2_end = floor((t_slant)/dt);

dphi_ref(dphi1_start:dphi1_end) = slant_ang/dphi_wind;
dphi_ref(dphi2_start:dphi2_end) = -slant_ang/dphi_wind;

for i=2:length(phi_ref)
    phi_ref(i) = phi_ref(i-1) + dphi_ref(i-1)*dt;
end


X_ref(1,:) = y_ref;
X_ref(2,:) = dy_ref;
X_ref(3,:) = phi_ref;
X_ref(4,:) = dphi_ref;

%% Controller
Q = eye(4);
R = 5.5;
Q(1,1) = 7.50;
Q(2,2) = 0.42;
Q(3,3) = 25;
Q(4,4) = 17.0;
[K,S,P] = lqr(A,B1,Q,R);

%% Simulation

X0 = [str_seg1_y,0,0,0]';
E = [0,0,0,0]';
traj_act = zeros(2, length(ts));
E_log = zeros(4,length(ts));
U_log = zeros(1,length(ts));
E_log(:,1) = [0,0,0,0]';


Edot = zeros(4,1);


for i=1:length(ts)
    curr_yaw = (phi_ref(i)+E(3));
    traj_act(1,i) = x_ref(i) - E(1)*sin(curr_yaw) - E(2)*sin(curr_yaw)*dt;
    traj_act(2,i) = y_ref(i) + E(1)*cos(curr_yaw) + E(2)*cos(curr_yaw)*dt;
    
    E_log(:,i) = E;
    U_log(i) = -K*E;
    Edot = A*E-B1*K*E+B2*dphi_ref(i);
    E = E + Edot*dt;
end


%% Display error analytics
fprintf('Max Abs. Lateral Error: %f m\n',max(abs(E_log(1,:))))
fprintf('Max Abs. Heading Error: %f rad\n',max(abs(E_log(3,:))))
fprintf('Max Steering Rate: %f rad/s\n',max(diff(U_log)/dt))

%% Plotting
figure;
plot(x_ref, y_ref)
hold on;
grid on;
plot(traj_act(1,:), traj_act(2,:));
title('Reference vs Tracking')
xlabel('X [m]')
ylabel('Y [m]')
legend('reference', 'actual')

figure;
plot(ts, E_log(1,:));
grid on;
title('Lateral Error')
xlabel('t [s]')
ylabel('y_{err} [rad]')
xline(t_strt1,'--','Lane Exit','LineWidth',1.50,'Color','b')
xline(t_strt1+1.0,'--','1 sec from Exit','LineWidth',1.50,'Color','b')
xline(t_slant,'--','Lane Entry','LineWidth',1.50,'Color','b')
xline(t_slant+1.0,'--','1 sec from Entry','LineWidth',1.50,'Color','b')
yline(0.002,'--','Settling Bound','LineWidth',1.5,'Color','g')
yline(-0.002,'--','Settling Bound','LineWidth',1.5,'Color','g')

figure;
plot(ts, U_log)
title('Steering Angle')
xlabel('t [s]')
ylabel('\delta_f [rad]')
grid on;

figure;
plot(ts(2:length(ts)),diff(U_log)/dt)
hold on;
grid on;
title('Steering Rate')
xlabel('t [s]')
ylabel('\delta_f'' [rad/s]')
xline(t_strt1,'--','Lane Exit','LineWidth',1.50,'Color','b')
xline(t_strt1+1.0,'--','1 sec from Exit','LineWidth',1.50,'Color','b')
xline(t_slant,'--','Lane Entry','LineWidth',1.50,'Color','b')
xline(t_slant+1.0,'--','1 sec from Entry','LineWidth',1.50,'Color','b')

figure;
plot(ts, E_log(3,:));
title('Heading Error')
xlabel('t [s]')
ylabel('\Phi_{err} [rad]')
xline(t_strt1,'--','Lane Exit','LineWidth',1.50,'Color','b')
xline(t_strt1+1.0,'--','1 sec from Exit','LineWidth',1.50,'Color','b')
xline(t_slant,'--','Lane Entry','LineWidth',1.50,'Color','b')
xline(t_slant+1.0,'--','1 sec from Entry','LineWidth',1.50,'Color','b')
yline(0.01,'--','Max Bound','LineWidth',1.5,'Color','r')
yline(-0.01,'--','Max Bound','LineWidth',1.5,'Color','r')
yline(0.0007,'--','Settling Bound','LineWidth',1.5,'Color','g')
yline(-0.0007,'--','Settling Bound','LineWidth',1.5,'Color','g')
grid on;


function U = control(X,B1,K,B2,dphi_des)
     U = -B1*K*X + B2*dphi_des;
end




