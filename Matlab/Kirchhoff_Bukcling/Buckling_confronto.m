clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% caratteristiche composito %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% moduli elastici della lamina e proprietà gemometriche
E_1 = 125e9; %(Pa)
E_2 = 12.5e9; %(Pa)
ni_12 = 0.38;
G_12 = 6.89e9; %(Pa)
t = 0.15e-3; %(m)


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
seq_theta1 = [0 90 0 90 0 90 0 90];
%seq_theta2 = [0 90 0 90 0 90 0 90];
seq_theta = [seq_theta1, fliplr(seq_theta1)];
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINIZIONE DEL PROBLEMA %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dati geometria 
a = 0.500; % [m] 
b = 0.150; % [m] 
t_lam = N*t; % [m] 
%Dati materiale 
%Dx = 9.27*10^4; % Rigidezza flessionale [N*mm] 
Dx= D(1,1);
Dy = D(2,2); % Rigidezza flessionale [N*m] 
Dxy = D(1,2); % Rigidezza flessionale [N*m] 
Ds = D(3,3); % Rigidezza flessionale [N*m] 
H = Dxy+2*Ds; 
n=1;
m_max=7;
N_x_cr_v=zeros(1,m_max);
for i=1:m_max
    N_x_cr_v(i)=(pi^2)*(Dx*((i/a)^2)+2*H*(n/b)^2+Dy*((n/b)^4)*((a/i)^2)); % [N/m]
end
[N_x_cr,m] = min(N_x_cr_v);
F_x_cr = N_x_cr_v*b;
n_passi = 700; 
X = linspace (0,a,n_passi); %N.B il numero di passi deve essere lo stesso 
Y = linspace (0,b,n_passi);
w = sin(m.*pi.*X./a).*sin(n.*pi.*Y./b)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%     PIASTRA ALLUMINIO             %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dati materiale 
E = 70e9; %% [MPa] 
ni = 0.33; % Modulo di Poisson 
k_c_v=zeros(1,m_max);
for i=1:m_max
    k_c_v(i)=(1/i^2)*((i/a)^2 + (n/b)^2)^2;
end
[k_c,m_al] = min(k_c_v);
D_al = (N_x_cr*m_al^2)/(pi*a*((m_al/a)^2+(n/b)^2))^2; % Rigidezza flesisonale [N*mm]
t_al = (12 * D_al *(1-ni^2)/E)^(1/3); % [mm]
n_passi = 700; 
X = linspace (0,a,n_passi); %N.B il numero di passi deve essere lo stesso 
Y = linspace (0,b,n_passi);
w_al = sin(m_al.*pi.*X./a).*sin(n.*pi.*Y./b)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Confronto tra i pesi   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vol_comp = a * b * t_lam ; % [m^3]
rho_comp = 1400; % [kg/m^3]
W_comp = Vol_comp * rho_comp;
Vol_al = a * b * t_al; % [m^3]
rho_al = 2800; % [kg/m^3]
W_al = Vol_al * rho_al;
rapp = W_comp/W_al;


disp (['Carico critico: ', num2str(N_x_cr),'N/m ']);
disp (['Spessore piastra in compodito: ', num2str(t_lam),'mm ']);
disp (['Spessore piastra in alluminio: ', num2str(t_al),'mm ']);
disp (['Peso piastra in alluminio: ', num2str(W_al),'kg']);
disp (['Peso pistra in composito: ', num2str(W_comp),'kg ']);
disp (['Rapporto pesi: ', num2str(rapp)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PLOT % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 
[X,Y] = meshgrid (X,Y); 
figure 
mesh(X,Y,w); 
xlabel([ '$' 'x\ (mm) ' '$'] ,'interpreter','latex'); 
ylabel([ '$' 'y\ (mm)' '$'] ,'interpreter','latex'); 
zlabel([ '$' 'w\ (mm)' '$'] ,'interpreter','latex'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PLOT % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 
figure 
mesh(X,Y,w_al); 
xlabel([ '$' 'x\ (mm) ' '$'] ,'interpreter','latex'); 
ylabel([ '$' 'y\ (mm)' '$'] ,'interpreter','latex'); 
zlabel([ '$' 'w\ (mm)' '$'] ,'interpreter','latex'); 