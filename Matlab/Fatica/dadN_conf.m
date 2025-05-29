clc;
clear;
close all;

% Constants and parameters
K_crit = 72.524;
A_k = 1;
B_k = 1;
C = 1.545e-10;
n = 3.284;
p = 0.5;
q = 1;
delta_K_0 = 3.187;
C_th = 1.5;
alpha = 72.524/36.262;
S_max_sigma_0 = 0.1;
a = 0.003;
a_0 = 0.0000;
R = -0.1;
K_max = linspace(7, K_crit-0.1);

% Calculating constants A_0, A_1, A_2, A_3
A_0 = (0.825 - 0.34*alpha + 0.05*alpha^2) * (cos(pi/2 * S_max_sigma_0))^(1/alpha);
A_1 = (0.415 - 0.071*alpha) * S_max_sigma_0;
A_3 = 2 * A_0 + A_1 - 1;
A_2 = 1 - A_0 - A_1 - A_3;

% Parameter lambda
lambda = 0.5;

% Conditional function for f based on R
if R >= 0
    f = max(R, A_0 + A_1*R + A_2*R^2 + A_3*R^3);
elseif -2 <= R && R < 0
    f = A_0 + A_1*R;
else
    f = A_0 - 2*A_1;
end

% Calculating delta_K_th
delta_K_th = delta_K_0 * sqrt(a/(a + a_0)) / (1 - f / ((1 - A_0) * (1 - R)))^(1 + C_th * R);

% Defining the delta_K range
delta_K = linspace(delta_K_th + 0.001, K_crit - 0.001);

% NASGRO equation
da_dN_nasgro = C * (((1 - f) / (1 - R)) * delta_K).^n .* (1 - delta_K_th ./ delta_K).^p ./ (1 - K_max / K_crit).^q;

dadN_Walker = [
    2	1.50491e-009
72.286	0.000196817
    ];

dadN_nas =[
3.08292	2.35577e-011
3.10368	1.11998e-010
3.13136	1.73113e-010
3.22132	3.14945e-010
3.42892	5.93046e-010
3.77493	1.10081e-009
4.46695	2.51187e-009
5.15896	4.64074e-009
6.54299	1.1854e-008
10.0031	5.76873e-008
16.9232	3.96695e-007
23.8434	1.44206e-006
30.7635	3.94908e-006
37.6837	9.32067e-006
44.6038	2.04073e-005
51.524	4.3911e-005
58.4441	0.000100005
65.3643	0.000289666
68.8244	0.000687023
71.7308	0.00491584
72.0768	0.0132792
72.2153	0.0397105
72.2222	0.044067
72.2775	0.352985
72.2845	2.79658
];



% Plotting NASGRO equation
figure(13);
loglog(delta_K, da_dN_nasgro, 'ob-');
axis([1 100 10^(-10) 10]);
xlabel('\DeltaK [Mpa*sqrt(m)]');
ylabel('da/dN [m/cycle]');
title('Confronto equazioni');
grid on;
hold on;
% Plotting NASGRO equation con AFGROW
loglog(dadN_nas(:,1), dadN_nas(:,2),'k-');

% Walker equation
da_dN_walker = C * delta_K.^n / (1 - R)^(n * (1 - lambda));
loglog(delta_K, da_dN_walker, 'or-');
% Walker equation con AFGROW
loglog(dadN_Walker(:,1), dadN_Walker(:,2),'k--','LineWidth',2);


% Paris equation
da_dN_paris = C * delta_K.^n;
loglog(delta_K, da_dN_paris, 'oy-');

% Legend and saving the figure
legend('NASGRO','NASGROWAFF', 'Walker','WalkerAFF', 'Paris');
hold off;
%saveas(gcf, 'figura5', 'jpg');
