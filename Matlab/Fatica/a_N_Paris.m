clc;
clear;
close all;

%% Materiale Al 2024T3
sigma_y = 365.4; % Mpa
E = 73.084e3; % Mpa
nu = 0.33;
G = E/(1 -2* nu); %Mpa
K_1c_ = 36.262; % Mpa*m^0.5
Ak = 1;   % coef corr spessore
Bk = 1;   % coef corr spessore
% PARAMETRI MODELLO DI WALKER
m = 3.284;
C = 1.545e-10;

%% Geometria del pannello e della cricca
W = 1; %m
t = 0.002; %m

% tenacità modificata per lo spessore
t_01 = 2.500 * (K_1c_/sigma_y)^2; %m

K_1c = K_1c_ * ( 1+Bk* exp(-Ak*(t/t_01)^2));

%% Esercizio 1 assegnato il valore della dimensione iniziale della cricca 
%% valuto la tenzione critica che fa propagare la cricca
a_in = 0.003; %m
Yi = 1 + 0.256*(a_in /W) - 1.152 * (a_in/W)^2 + 12.2*(a_in/W)^3;
sigma_1cr = K_1c/(Yi*((pi*a_in)^0.5));

disp(['Il carico statico che porta alla propagazione della cricca è: ',num2str(sigma_1cr),'MPa'])
%disp(['Poichè è più grande del carico di snervamento, la cricca non propaga staticamente'])
%% Esercizio 2 assegnato il valore della tenzione massima 
%% valuto la dimensione della cricca a cui inzia la propagazione instabile
sigma_max = 80;
a_cr = (1/pi)*(K_1c/(sigma_max))^2;
disp(['Senza tener conto del Correction Factor Y'])
disp(['La dimensione della cricca a partire dalla quale c''è propagazione: ',num2str(a_cr),'m'])
Y_f_errato = 1 + 0.256*(a_cr /W) - 1.152 * (a_cr/W)^2 + 12.2*(a_cr/W)^3;
disp(['Valore del fattore correttivo: ',num2str(Y_f_errato)])
%% Metodo iterativo per trovare Y
err = 2;
tol = 1e-6;
i = 1;
Y_f_cor = 1;

while err > tol
a_cr_cor = (1/pi)*(K_1c/(Y_f_cor*sigma_max))^2;
Y_f_new = 1 + 0.256*(a_cr_cor /W) - 1.152 * (a_cr_cor/W)^2 + 12.2*(a_cr_cor/W)^3;
Y_f_cor = (Y_f_new + Y_f_cor)/2;
err =abs(Y_f_cor - Y_f_new);
i = i+1;
end

trova_a = @(a) K_1c - (sigma_max * ...
    (1 + 0.256*(a_cr_cor /W) - 1.152 * (a_cr_cor/W)^2 + 12.2*(a_cr_cor/W)^3) *...
    ((pi*a)^0.5));
a_tent = fzero(trova_a,0.19);
disp(['Considerando il Correction Factor Y'])
disp(['La dimensione della cricca a partire dalla quale c''è propagazione: ',num2str(a_cr_cor),'m'])
disp(['Valore del fattore correttivo: ',num2str(Y_f_cor)])
%% Esercizio 3 assegnato il valore della tenzione massima 
%% valutata la dimensione della cricca a cui inzia la propagazione instabile
%% calcolo numero di cicli necessari per raggiungere quella dimensione della cricca
%% calcolo come varia la dimensione della cricca al crescere del numero di cicli

a_fin = a_cr;

Nf = (a_fin^(1-(m/2)) - a_in^(1-(m/2)))/(C*(1-(m/2))*(sigma_max^m)*(pi^(m/2)));

N = 1 : 250: Nf;
a_N = (a_in^(1-(m/2))+ (C*(1-(m/2))*(sigma_max^m)*(pi^(m/2))*N)).^(1/(1-(m/2)));
a_N_prov = (a_in^(1-(m/2))+ (C*(1-(m/2))*(sigma_max^m)*(pi^(m/2))*33700)).^(1/(1-(m/2)));

dadN = @(N,a) C*(sigma_max*...
    (1 + 0.256*(a_cr /W) - 1.152 * (a_cr/W)^2 + 12.2*(a_cr/W)^3)*...
    ((pi*a)^0.5))^m;


Nff = 2.5e4;
options = odeset('Events', @myevent);  % Use the correct event function handle

% Solve the ODE with event function
[N_v, a_v, N_e, a_e, i_e] = ode45(dadN, [0 Nff], a_in, options);


%% AFGROW

Cycles_C = [...
0	0.003
17200	0.00817642
22400	0.0134206
25100	0.0187867
26800	0.0242724
27900	0.0293952
28800	0.0350666
29500	0.0408507
30000	0.0460232
30500	0.0523837
30900	0.0586068
31200	0.0641343
31500	0.0706049
31800	0.0782593
32000	0.0841917
32200	0.0909435
32400	0.0986951
32600	0.107687
32800	0.118248
32898	0.12417
32994	0.13054
33090	0.137575
33185	0.145311
33281	0.15406
33376	0.163842
33471	0.175007
33566	0.187926
33600	0.193214
33660	0.202981
33700	0.210435
];



%% grafica


figure(2)
    plot(N,a_N,"Color",[0 0.4470 0.7410],'LineWidth',1.5);
    hold on,
    grid on,
    plot(Cycles_C(:,1),Cycles_C(:,2),"Color",[0.8500 0.3250 0.0980],'LineWidth',1.5);
    plot(N_v,a_v,"Color",[0.4940 0.1840 0.5560],'LineWidth',1.5);
    title('Modello di Paris');
ylabel('$a$','Interpreter','latex','FontSize',24);
xlabel('N ','Interpreter','latex','FontSize',24);
lgd = legend('Affgrow',' Soluzione Numerica', 'SolconBcorrection');
%lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;




function [position, isterminal, direction] = myevent(N, a)
    position = 0.2091 - a; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end

% function [position,isterminal,direction] = myevent(N,a)
%   sigma = 80;
%   Kcr = 72.286;
%   W=1;
%   Y = (1 + 0.256*(a /W) - 1.152 * (a/W)^2 + 12.2*(a/W)^3);
%   K = sigma * Y * (a * pi);
%   %position = Kcr-K; % The value that we want to be zero
%   position = 0.2599-a
%   isterminal = 0;  % Halt integration 
%   direction = 0;   % The zero can be approached from either direction
% end


