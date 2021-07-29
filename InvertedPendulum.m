%% Control Systems Design
% Dimitris Dimos [031 17 165]
% <dimitris.dimos647@gmail.com>
% 6th Semester - Flow S
% Lab Exercise - Inverted Pendulum


%% -------------  Question A  -------------

% Matrices
A = [0 1 0 0; 20.6 0 0 0; 0 0 0 1; -.5 0 0 0];
B = [0; -1; 0; .5];
C = [1 0 0 0; 0 0 1 0];
D = [0; 0];

x0 = [-.2; -.06; .01; .3];  % initial conditions

% State-space model
states = {'theta' 'theta_dot' 'x' 'x_dot'};
inputs = {'u'};
outputs = {'theta'; 'x'};
system = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

% System Controllability
control_matrix = ctrb(system);
controllability = rank(control_matrix);

% Design requirements
Ts = 1.8;      % settling time (<= 2 seconds)
zeta = 0.5;  % damping ratio

wn = 4/(Ts*zeta);   % natural frequency
wd = wn*sqrt(1-zeta^2);

% Desired poles locations
pole1 = -25;
pole2 = -26;
pole3 = -zeta*wn + wd*i;
pole4 = -zeta*wn - wd*i;
p = [pole1 pole2 pole3 pole4];

K = place(A, B, p);   % K matrix to move poles


%% -------------  Question B  -------------

% Riccati Equation Solution
[K_lqr, P, e] = lqr(system, eye(4), 1, 0); % solution P matrix and optimal K_lqr gain

%% -------------  Question C  -------------

% We solve for u(t): 0 = (A-BK)x(tf) + Bu => -Bu = (A-BK)x(tf)
x_tf = [0; 0; 1; 0];
Bu_tf = -(A - B*K)*x_tf;
utf = B\Bu_tf;  % we find that utf = K(3)

%% -------------  Question D  -------------

observability_matrix = obsv(system);
observability = rank(observability_matrix);
obv_poles = [-10 -11 -12 -13];
L = place(A', C', obv_poles)';

%% -------------  Question E  -------------

A_last = [0 1 0 0; 20.9 0 0 0; 0 0 0 1; -.8 0 0 0];
system_last = ss((A_last-B*K),[0;0;0;0],C,D,'statename',states,'inputname',inputs,'outputname',outputs);
