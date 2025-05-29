clc;
clear;
close all;

%% moduli elastici della lamina e proprietà gemometriche
E_1 = 125e9; %(Pa)
E_2 = 12.5e9; %(Pa)
ni_12 = 0.38;
G_12 = 6.89e9; %(Pa)
t = 0.15e-3; %(m)
L= 0.5; %(m)


%% matrice Q
Q_11 = E_1/(1-(E_2/E_1)*(ni_12)^2);
Q_12 = ni_12*E_2/(1-(E_2/E_1)*(ni_12)^2);
Q_22 = E_2/(1-(E_2/E_1)*(ni_12)^2);
Q_66 = G_12;

Q=[Q_11 Q_12 0;
    Q_12 Q_22 0;
    0 0 Q_66];

%% calcolo matrici di rotazione

T_sigma = @(theta) [(cos(theta))^2  (sin(theta))^2  -2*cos(theta)*sin(theta);
    (sin(theta))^2  (cos(theta))^2  2*cos(theta)*sin(theta);
    cos(theta)*sin(theta)  -cos(theta)*sin(theta)  (cos(theta))^2-(sin(theta))^2];


T_eps = @(theta) [(cos(theta))^2  (sin(theta))^2  -cos(theta)*sin(theta);
    (sin(theta))^2  (cos(theta))^2  cos(theta)*sin(theta);
    2*cos(theta)*sin(theta)  -2*cos(theta)*sin(theta)  (cos(theta))^2-(sin(theta))^2];


%% calcolo matrici Q nel riferimento globale

Q_glob = @(theta) (T_sigma(theta) * Q) / T_eps(theta);

%% sequenza di laminazione 
seq_theta1 = [90 45 -45 0 90 45 -45 0 90 45 -45 0];
seq_theta = [seq_theta1, -fliplr(seq_theta1)];%fliplr(seq_theta1)];
seq_theta_rad = convang(seq_theta,'deg','rad');

N = length(seq_theta_rad);
z_vec = t * linspace((-N/2),(N/2),N+1);
z_mean= z_vec(2:end)-t/2;
%% costruzione matrici A, B, D
A = zeros(3);
B = zeros(3);
D = zeros(3);

for i=1:N
    A = A + Q_glob(seq_theta_rad(i)) * (z_vec(i+1) - z_vec(i));
    B = B + Q_glob(seq_theta_rad(i)) * ((z_vec(i+1))^2 - (z_vec(i))^2)/2;
    D = D + Q_glob(seq_theta_rad(i)) * ((z_vec(i+1))^3 - (z_vec(i))^3)/3;
end
ABBD = [A B; B D];