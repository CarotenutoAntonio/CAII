clc;
clear;
close all;

%% moduli elastici della lamina e proprietà gemometriche
E_1 = 125e9; %(Pa)
E_2 = 12.5e9; %(Pa)
ni_12 = 0.38;
G_12 = 6.89e9; %(Pa)
t = 0.15e-3; %(m)

%% carichi esterni apllicati 
Delta_P = 0.18e6; %(Pa) carico di pressurizzazione
M_T = 4.5e6; %(Nm)
R = 2; %(m) raggio di fusoliera
M_S = 0.5; %(m) Margine di sicurezza assegnato
fat_ampl = 1+M_S; % per il calcolo
%valuto carichi per unità di lunghezza
M_x = 0; %(Nm)
M_y = 0; %(Nm)
M_xy = 0;%(Nm)
N_x = Delta_P * R/2 ; %(N/m)
N_y = Delta_P * R; %(N/m)
N_xy = M_T/(2*pi*R^2); %(N/m)

N_M_ul = [N_x N_y N_xy M_x M_y M_xy];


%% sequenza di laminazione 
%% scelgo un laminato simmetrico e trasversalmente isotropo
seq_theta1 = [0 45 -45 90];

% Definisco carichi ultimi laminato usando carpet plot nota la sequenza
sigma_xLT = 172e6; %(Pa)
sigma_yLT = 172e6; %(Pa)
tauxy_L = 96e6; %(Pa)

%% calcolo spessore e quindi numero di lamine
N_xt = ceil(N_x*fat_ampl/(sigma_xLT*8*t));
N_yt = ceil(N_y*fat_ampl/(sigma_yLT*8*t));
N_st = ceil(N_xy*fat_ampl/(tauxy_L*8*t));
M = max([N_xt, N_yt, N_st]);
seq_theta2=seq_theta1;
for i=1:M-1
    % if i==1
    % seq_theta2=[seq_theta2, fliplr(seq_theta1)];
    % else
        seq_theta2=[seq_theta2, seq_theta1];
    %end

end
seq_theta = [ 45, -45 seq_theta2, fliplr(seq_theta2), 45, -45];
%seq_theta = [seq_theta2, fliplr(seq_theta2)];
seq_theta_rad = convang(seq_theta,'deg','rad');
N = length(seq_theta_rad);
t_lam= N*t;
lamine = 1:1:N;
fprintf('Il numero di lamine è %f \n',N);
fprintf('Lo spessore è %f m \n ',t_lam);
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

%% costruzione matrici A, B, D
z_vec = t * linspace((-N/2),(N/2),N+1);
z_mean= z_vec(2:end)-t/2;
A = zeros(3);
B = zeros(3);
D = zeros(3);

for i=1:N
    A = A + Q_glob(seq_theta_rad(i)) * (z_vec(i+1) - z_vec(i));
    B = B + Q_glob(seq_theta_rad(i)) * ((z_vec(i+1))^2 - (z_vec(i))^2)/2;
    D = D + Q_glob(seq_theta_rad(i)) * ((z_vec(i+1))^3 - (z_vec(i))^3)/3;
end
ABBD = [A B; B D];
ABBD_inv = inv(ABBD);

%% moduli elastici effettivi
E_x = 1/(ABBD_inv(1,1)*t*N);
E_y = 1/(ABBD_inv(2,2)*t*N);
G_xy = 1/(ABBD_inv(3,3)*t*N);
ni_xy = - ABBD_inv(1,2) / ABBD_inv(1,1);

%% moduli di bending
E_bx = 12/(ABBD_inv(4,4)*(t*N)^3);
E_by = 12/(ABBD_inv(5,5)*(t*N)^3);


%% calcolo deformazioni
eps_k = ABBD\N_M_ul(:);
eps_mat_glob = zeros(3,N);
sigma_mat_glob = zeros(3,N);
eps_mat_loc = zeros(3,N);
sigma_mat_loc = zeros(3,N);

for i=1:N

    eps_mat_glob(:,i) = eps_k(1:3) + z_mean(i)*eps_k(4:6);
    sigma_mat_glob(:,i) = Q_glob(seq_theta_rad(i))*eps_mat_glob(:,i);

    eps_mat_loc(:,i) = T_eps(seq_theta_rad(i))\eps_mat_glob(:,i);
    sigma_mat_loc(:,i) = T_sigma(seq_theta_rad(i))\sigma_mat_glob(:,i);

end

%% Indici di Failure
%% Definisco carichi ultimi lamina usando carpet plot
%% calcolo margini di stabilità

sigma_1LT = 434e6; %(Pa)
sigma_2LT = 41.4e6; %(Pa)
sigma_1LC = -331e6; %(Pa)
sigma_2LC = -27.6e6; %(Pa)
tau_L = 34.5e6; %(Pa)

eps_1LT =  sigma_1LT / E_1;
eps_2LT =  sigma_2LT / E_2;
eps_1LC =  sigma_1LC / E_1;
eps_2LC =  sigma_2LC / E_2;
gam_L = tau_L / G_12;

Fail_ind= zeros(7,1);
I = zeros(6,1);

MaxStress1 = zeros(1,N);
MaxStress2 = zeros(1,N);
MaxStressS = zeros(1,N);

MaxStrain1 = zeros(1,N);
MaxStrain2 = zeros(1,N);
MaxStrainS = zeros(1,N);
TH = zeros(1,N);

%% max stress e max strain
for i=1:N

if sigma_mat_loc(1,i)<0
    MaxStress1(i) = sigma_1LC/sigma_mat_loc(1,i);
    % costruisco un termine per applicare Tsai-Hill
    TH(1,i) = (sigma_mat_loc(1,i)*sigma_mat_loc(2,i))/(sigma_1LC^2);
else 
    MaxStress1(i) = sigma_1LT/sigma_mat_loc(1,i);
     % costruisco un termine per applicare Tsai-Hill
    TH(1,i) = (sigma_mat_loc(1,i)*sigma_mat_loc(2,i))/(sigma_1LT^2);
end


if sigma_mat_loc(2,i)<0
    MaxStress2(i) = sigma_2LC/sigma_mat_loc(2,i);
else 
    MaxStress2(i) = sigma_2LT/sigma_mat_loc(2,i);
end

if eps_mat_loc(1,i)<0
    MaxStrain1(i) = eps_1LC/eps_mat_loc(1,i);
else 
    MaxStrain1(i) = eps_1LT/eps_mat_loc(1,i);
end

if eps_mat_loc(2,i)<0
    MaxStrain2(i) = eps_2LC/eps_mat_loc(2,i);
else 
    MaxStrain2(i) = eps_2LT/eps_mat_loc(2,i);
end


end

MaxStressS = abs(tau_L./sigma_mat_loc(3,:));
MaxStrainS = abs(gam_L./eps_mat_loc(3,:));

%% Max stress

Maxstres=min([MaxStress1',MaxStress2',MaxStressS'],[],2);

%% Max strain

Maxstrain=min([MaxStrain1',MaxStrain2',MaxStrainS'],[],2);


%% Tsai Hill
TsaiHill = (1./(((1./MaxStress1).^2)+((1./MaxStress2).^2)+((1./MaxStressS).^2) - TH)).^(1/2);

%% Tsai-Wu
f1 = (1/(sigma_1LT)) + (1/(sigma_1LC));
f11 = -1/(sigma_1LT*sigma_1LC);
f2 = (1/(sigma_2LT)) + (1/(sigma_2LC));
f22 = -1/(sigma_2LT*sigma_2LC);
f66 = 1/(tau_L^2);
f12 = -(1/2)*((f11*f22)^(1/2));

BB =  f1*sigma_mat_loc(1,:) + f2*sigma_mat_loc(2,:);
AA =  f11*(sigma_mat_loc(1,:).^2) + f22*(sigma_mat_loc(2,:).^2) + ...
    f66*(sigma_mat_loc(3,:).^2) + 2*f12*(sigma_mat_loc(1,:).*(sigma_mat_loc(2,:)));
C = -1;

TsaiWu = (-BB + ((BB.^2 - 4*AA.*C)).^(1/2))./(2*AA);



%% calcolo indice di sicurezza minimo
[Fail_ind(1), I(1)] = min(MaxStress1);
[Fail_ind(2), I(2)] = min(MaxStress2);
[Fail_ind(3), I(3)] = min(MaxStressS);
[Fail_ind(4), I(4)] = min(MaxStrain1);
[Fail_ind(5), I(5)] = min(MaxStrain2);
[Fail_ind(6), I(6)] = min(MaxStrainS);
[Fail_ind(7), I(7)] = min(TsaiHill);
[Fail_ind(8), I(8)] = min(TsaiWu);

%Calcolo Failure index minimo e mostro
[Fail_index,pos]= min (Fail_ind);
fprintf('Indice di failure è %f \n',Fail_index);
strin = ["MaxStress1", "MaxStress2", "MaxStressS", ...
         "MaxStrain1", "MaxStrain2", "MaxStrainS", "TsaiHill", "TsaiWu"];

d=strin(pos);
disp('Si è ottenuto con il criterio di ');
disp(d);
fprintf('nella lamina  %f \n',I(pos));

%% failure index tabella
tab=table(lamine',seq_theta',Maxstres,Maxstrain,TsaiHill',TsaiWu');
tab.Properties.VariableNames=[{'lamina'},{'seq_theta'},{'MaxStress'},{'MaxStrain'},{'TsaiHill'},{'TsaiWu'}];
% Calculate the size of the table based on content
tableWidth = 120 * size(tab, 2); % Width per column
tableHeight = 22 * size(tab, 1); % Height per row
% Creating a figure with a uitable
figure;
uitable('Data', tab{:,:}, 'ColumnName', tab.Properties.VariableNames, ...
'RowName', [], 'Position', [100, 100, tableWidth, tableHeight]);
sgtitle('Indici di Failure, con diversi criteri');

%% stress e deformazioni globali tabella
% Creazione della tabella
tab = table(lamine',seq_theta', sigma_mat_glob(1,:)', sigma_mat_glob(2,:)', sigma_mat_glob(3,:)', ...
    eps_mat_glob(1,:)', eps_mat_glob(2,:)', eps_mat_glob(3,:)');

% Impostazione dei nomi delle colonne con simbolo greco σ
tab.Properties.VariableNames = [{'lamina'},{'seq_theta'},{'σx'}, {'σy'}, {'τxy'}, {'εx'}, {'εy'}, {'γxy'}];

% Calcolare la larghezza e l'altezza della tabella in base al contenuto
tableWidth = 120 * size(tab, 2);  % Larghezza per colonna
tableHeight = 22 * size(tab, 1);  % Altezza per riga

% Creare una figura con una uitable
figure;

% Creazione della tabella
uitable('Data', tab{:,:}, 'ColumnName', tab.Properties.VariableNames, ...
    'RowName', [], 'Position', [100, 100, tableWidth, tableHeight]);

% Titolo con simbolo LaTeX
sgtitle('Stress e deformazioni nel riferimento globale', 'Interpreter', 'latex');


%% stress e deformazioni locali tabella
% Creazione della tabella
tab = table(lamine',seq_theta', sigma_mat_loc(1,:)', sigma_mat_loc(2,:)', sigma_mat_loc(3,:)', ...
    eps_mat_loc(1,:)', eps_mat_loc(2,:)', eps_mat_loc(3,:)');

% Impostazione dei nomi delle colonne con simbolo greco σ
tab.Properties.VariableNames = [{'lamina'},{'seq_theta'},{'σ1'}, {'σ2'}, {'τ12'}, {'ε1'}, {'ε2'}, {'γ12'}];

% Calcolare la larghezza e l'altezza della tabella in base al contenuto
tableWidth = 120 * size(tab, 2);  % Larghezza per colonna
tableHeight = 22 * size(tab, 1);  % Altezza per riga

% Creare una figura con una uitable
figure;

% Creazione della tabella
uitable('Data', tab{:,:}, 'ColumnName', tab.Properties.VariableNames, ...
    'RowName', [], 'Position', [100, 100, tableWidth, tableHeight]);

% Titolo con simbolo LaTeX
sgtitle('Stress e deformazioni nel riferimento locale', 'Interpreter', 'latex');


