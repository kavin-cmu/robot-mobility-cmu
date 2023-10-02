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
Vx = 100;

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

% Measurement Matrix
C = [1 0 0 0;
     0 0 1 0];

%% Reference Path Generation
strt_seg1_dur = 1.0;
strt_seg1_len = Vx*strt_seg1_dur;

strt_seg2_dur = 1.0;
strt_seg2_len = Vx*strt_seg2_dur;

circ_seg1_rad = 1000;
circ_seg1_dur = 5.0;
circ_seg1_dir = 1;

circ_seg2_rad = 500;
circ_seg2_dur = 5.0;
circ_seg2_dir = -1;

init_pos = [0,0]';
init_phi = 0.0;

%% Time Parametrization
dt = 0.01;
tf = strt_seg1_dur+strt_seg2_dur+circ_seg1_dur+circ_seg2_dur;
t_strt1 = strt_seg1_dur;
t_circ1 = t_strt1 + circ_seg1_dur;
t_strt2 = t_circ1 + strt_seg2_dur;
t_circ2 = t_strt2 + circ_seg2_dur;
ts1 = linspace(0,t_strt1,floor(t_strt1/dt));
ts2 = linspace(t_strt1+dt,t_circ1,floor((circ_seg1_dur)/dt));
ts3 = linspace(t_circ1+dt,t_strt2,floor((strt_seg2_dur)/dt));
ts4 = linspace(t_strt2+dt,t_circ2,floor((circ_seg2_dur)/dt));

ts = [ts1, ts2, ts3, ts4];

dy_ref = zeros(1,length(ts));
phi_ref = zeros(1,length(ts));
dphi_ref = zeros(1,length(ts));
X_ref = zeros(4,length(ts));

path_ref = zeros(2,length(ts));

%straight segment 1
path_ref(:,1:length(ts1)) = init_pos + [Vx*cos(init_phi), Vx*sin(init_phi)]'*ts1;
dy_ref(1:length(ts1)) = Vx*sin(init_phi);

%Circular segment 1
circ1_end_phi = init_phi+Vx*circ_seg1_dur/circ_seg1_rad;
theta1 = init_phi+Vx*linspace(0, circ_seg1_dur, length(ts2))/circ_seg1_rad;
circ1_pts = circ_seg1_rad*[cos(theta1-pi/2); sin(theta1-pi/2)];
circ1_trans = path_ref(:,length(ts1)) - circ1_pts(:,1);
circ1_pts = circ1_pts + circ1_trans + [Vx*cos(init_phi), Vx*sin(init_phi)]'*dt;
path_ref(:,length(ts1)+1:length(ts1)+length(ts2)) = circ1_pts;
dy_ref(length(ts1)+1:length(ts1)+length(ts2)) = Vx*sin(theta1);

%Straight Segment 2
path_ref(:,length(ts1)+length(ts2)+1:length(ts1)+length(ts2)+length(ts3)) = path_ref(:,length(ts1)+length(ts2)) + [Vx*cos(circ1_end_phi), Vx*sin(circ1_end_phi)]'*(ts3-t_circ1);
dy_ref(length(ts1)+length(ts2)+1:length(ts1)+length(ts2)+length(ts3)) = Vx*sin(circ1_end_phi);

%Circular segment 2
circ2_end_phi = circ1_end_phi-Vx*circ_seg2_dur/circ_seg2_rad;
theta2 = circ1_end_phi - Vx*linspace(0, circ_seg2_dur, length(ts4))/circ_seg2_rad;
circ2_pts = circ_seg2_rad*[cos(theta2+pi/2); sin(theta2+pi/2)];
circ2_trans = path_ref(:,length(ts3)+length(ts2)+length(ts1))+[Vx*cos(circ1_end_phi), Vx*sin(circ1_end_phi)]'*dt - circ2_pts(:,1);
circ2_pts = circ2_pts + circ2_trans;
path_ref(:,length(ts1)+length(ts2)+length(ts3)+1:length(ts1)+length(ts2)+length(ts3)+length(ts4)) = circ2_pts;
dy_ref(length(ts1)+length(ts2)+length(ts3)+1:length(ts1)+length(ts2)+length(ts3)+length(ts4)) = Vx*sin(theta2);

for i=1:length(phi_ref)-1
    phi_ref(i) = atan2(path_ref(2,i+1)-path_ref(2,i),path_ref(1,i+1)-path_ref(1,i));
end

phi_ref(length(phi_ref)) = phi_ref(length(phi_ref)-1);

for i=1:length(dphi_ref)-1
    dphi_ref(i) = (phi_ref(i+1)-phi_ref(i))/dt;
end

dphi_ref(length(dphi_ref)-1) = dphi_ref(length(dphi_ref)-2);
dphi_ref(length(dphi_ref)) = dphi_ref(length(dphi_ref)-2);

X_ref(1,:) = path_ref(2,:);
X_ref(2,:) = dy_ref;
X_ref(3,:) = phi_ref;
X_ref(4,:) = dphi_ref;

%% Controller
Q = eye(4);
R = 1.0;
Q(1,1) = 15;
Q(2,2) = 0.005;
Q(3,3) = 30;
Q(4,4) = 0.005;
[K,S,P] = lqr(A,B1,Q,R);

%% Simulation 1

X0 = [init_pos(1),init_pos(2),init_phi,0]';
X = [0,0,0,0]';
E = [0,0,0,0]';
X_log = zeros(4,length(ts));
E_log = zeros(4,length(ts));
U_log = zeros(1,length(ts));
Edot = zeros(4,1);
traj_act = zeros(2, length(ts));

for i=1:length(ts)
    X_log(:,i) = X_ref(:,i)+E;
    curr_yaw = (phi_ref(i)+E(3));
    traj_act(1,i) = path_ref(1,i) - E(1)*sin(curr_yaw) - E(2)*sin(curr_yaw)*dt;
    traj_act(2,i) = path_ref(2,i) + E(1)*cos(curr_yaw) + E(2)*cos(curr_yaw)*dt;
    E_log(:,i) = E;
    U_log(i) = -K*E;
    Edot = A*E-B1*K*E+B2*dphi_ref(i);
    E = E + Edot*dt;
end

% Print analysed variables
fprintf('Max Abs. Lateral Error: %f m\n',max(abs(E_log(1,:))))
fprintf('Max Abs. Heading Error: %f rad\n',max(abs(E_log(3,:))))
fprintf('Max Steering Rate: %f rad/s\n',max(diff(U_log)/dt))

%% Plotting

figure;
plot(ts, abs(E_log(1,:)));
title('Lateral Error')
xlabel('t [s]')
ylabel('y_{err} [m]')
yline(0.01,'--','Bound','LineWidth',2.0,'Color','r')
grid on;

figure;
plot(ts, U_log)
title('Steering Angle')
xlabel('t [s]')
ylabel('\delta_f [rad]')
grid on;

figure;
plot(ts(1:length(ts)-1),diff(U_log)/dt)
title('Steering Rate')
xlabel('t [s]')
ylabel('\delta_f'' [rad/s]')
grid on;
grid on;

figure
plot(ts, abs(E_log(3,:)));
grid on;
title('Heading Error')
xlabel('t [s]')
ylabel('\Phi_{err} [rad]')
yline(0.01,'--','Bound','LineWidth',2.0,'Color','r')

figure;
plot(path_ref(1,:),path_ref(2,:));
hold on;
plot(path_ref(1,:), X_log(1,:))
grid on;
title('Reference vs Actual Path')
xlabel('x [m]')
ylabel('y [m]')
legend('Reference', 'Actual')


function U = control(X,B1,K,B2,dphi_des)
     U = -B1*K*X + B2*dphi_des;
end




