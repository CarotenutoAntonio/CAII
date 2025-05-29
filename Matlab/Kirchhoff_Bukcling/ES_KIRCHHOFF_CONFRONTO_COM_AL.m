clc; clear all; close all; 
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
seq_theta1 = [90 90 90 90 90 90 90 90];
%seq_theta2 = [0 90 0 90 0 90 0 90];
seq_theta = [seq_theta1, seq_theta1];%fliplr(seq_theta1)];
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
a = 400; % [mm] 
b = 400; % [mm] 
t_ply = t*16*10^3; % [mm] 
%Dati materiale 
%Dx = 9.27*10^4; % Rigidezza flessionale [N*mm] 
Dx= D(1,1)*10^3;
Dy = D(2,2)*10^3; % Rigidezza flessionale [N*mm] 
Dxy = D(1,2)*10^3; % Rigidezza flessionale [N*mm] 
Ds = D(3,3)*10^3; % Rigidezza flessionale [N*mm] 
H = Dxy+2*Ds; 
%Carico costante 
q = -0.0001; % carico di pressione [MPa] 
n_passi = 700; 
X = linspace (0,a,n_passi); %N.B il numero di passi deve essere lo stesso 
Y = linspace (0,b,n_passi); 
% Inizializzazione 
w_tot = zeros(n_passi); 
w1 = zeros(n_passi); 
n_interazioni = 13; 
n_m = ceil(n_interazioni/2); 
vettore_w = zeros(1,n_m); 
h = 1; 
M = linspace(1,n_m,n_m); 

% Doppia serie 
for n=1:2:n_interazioni 
for m = 1:2:n_interazioni 
w = (16*q)/(pi^6)*((1/(m*n))*sin(m.*pi.*X./a).*sin(n.*pi.*Y./b)')/... 
(Dx*(m/a)^4+2*H*(m/a)^2*(n/b)^2+Dy*(n/b)^4); 
w1 = w;
w_tot = w_tot+w; 
if (m == n)
minimo = min(min(w_tot)); 
vettore_w(h)=minimo; 
h = h+1; 
disp (['w', num2str(m), '(a/2,b/2) = ', num2str(minimo)]); 
end 
end 
end

disp (['Massimo spostamento in entrambi casi: ', num2str(minimo),'mm ']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PLOT % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 
[X,Y] = meshgrid (X,Y); 
figure 
mesh(X,Y,w_tot); 
title('Deformata piastra in composito')
xlabel([ '$' 'x\ (mm) ' '$'] ,'interpreter','latex'); 
ylabel([ '$' 'y\ (mm)' '$'] ,'interpreter','latex'); 
zlabel([ '$' 'w\ (mm)' '$'] ,'interpreter','latex'); 
figure 
plot (M,vettore_w,'o-b', 'Linewidth' ,1.5); 
grid on; 
hold on 
xlabel([ '$' 'm\ ' '$'] ,'interpreter','latex'); 
ylabel([ '$' 'w(x,y)\ (mm)' '$'] ,'interpreter','latex');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%   PIASTRA ALLUMINIO            %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
 %% Dati geometria 
a = 400; % [mm] 
b = 400; % [mm] 
%t_tot = 2.4; % [mm] 
 
%% Dati materiale 
E = 70000; %% [MPa] 
ni = 0.33; % Modulo di Poisson 
coef_analitico = 0.04343;
t_tot = (abs(coef_analitico*q*(a^4)/(E*minimo)))^(1/3); % [mm]
% t_tot = 2.4;
D = (E*t_tot^3)/(12*(1-ni^2)); % Rigidezza flesisonale [N*mm]


%% Discretizzazione 
n_passi = 700; 
X = linspace (0,a,n_passi); %il numero di passi deve essere lo stesso 
Y = linspace (0,b,n_passi); 

%% Inizializzazione 
w_tot = zeros(n_passi); 
w1 = zeros(n_passi); 
n_interazioni = 13; 
n_m = ceil(n_interazioni/2); 
vettore_w = zeros(1,n_m);


h = 1; 
M = linspace(1,n_m,n_m); 

%%Doppia serie 
for n=1:2:n_interazioni 
for m = 1:2:n_interazioni 
w = (16*q)/(D*pi^6)*((1/(m*n))*sin(m.*pi.*X./a).*sin(n.*pi.*Y./b)')/... 
 (((m/a)^2+(n/b)^2)^2); 
w1 = w; 
w_tot = w_tot+w; 
if (m == n) 
minimo = min(min(w_tot)); 
vettore_w(h)=minimo; 
h = h+1; 
disp (['w', num2str(m), '(a/2,b/2) = ', num2str(minimo)]); 
end 
end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% PLOT % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 
[X,Y] = meshgrid (X,Y); 
figure 
mesh(X,Y,w_tot); 
title('Deformata piastra in alluminio')
xlabel([ '$' 'x\ (mm) ' '$'] ,'interpreter','latex'); 
ylabel([ '$' 'y\ (mm)' '$'] ,'interpreter','latex'); 
zlabel([ '$' 'w\ (mm)' '$'] ,'interpreter','latex'); 
figure 
plot (M,vettore_w,'o-b', 'Linewidth' ,1.5); 
grid on; 
hold on 
xlabel([ '$' 'm\ ' '$'] ,'interpreter','latex'); 
ylabel([ '$' 'w(x,y)\ (mm)' '$'] ,'interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Confronto tra i pesi   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vol_comp = a * b * t_ply * 10^-9; % [m^3]
rho_comp = 1400; % [kg/m^3]
W_comp = Vol_comp * rho_comp;

Vol_al = a * b * t_tot * 10^-9; % [m^3]
rho_al = 2800; % [kg/m^3]
W_al = Vol_al * rho_al;

rapp = W_comp/W_al;
disp (['Massimo spostamento in entrambi casi: ', num2str(minimo),'mm ']);
disp (['Spessore piastra in compodito: ', num2str(t_ply),'mm ']);
disp (['Spessore piastra in alluminio: ', num2str(t_tot),'mm ']);
disp (['Peso piastra in alluminio: ', num2str(W_al),'mm ']);
disp (['Peso pistra in composito: ', num2str(W_comp),'kg ']);
disp (['Rapporto pesi: ', num2str(rapp)]);