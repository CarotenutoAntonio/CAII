clc;
clear;
close all;

E_1 = 125e9; %(Pa)
E_2 = 12.5e9; %(Pa)
nu_12 = 0.38;
G_12 = 6.89e9; %(Pa)
t = 0.15e-3; %(m)
L= 0.5; %(m)


%% matrice Q
Q_11 = E_1/(1-(E_2/E_1)*(nu_12)^2);
Q_12 = nu_12*E_2/(1-(E_2/E_1)*(nu_12)^2);
Q_22 = E_2/(1-(E_2/E_1)*(nu_12)^2);
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



%% calcolo matrice Q nel riferimento globale e matrice S

Q_glob = @(theta) (T_sigma(theta) * Q) / T_eps(theta);

S_glob = @(theta) inv(Q_glob(theta));

%% plot
N = 20;
theta_vec = linspace(0,pi/2,N);
E_x_vec = zeros(1,N);
E_y_vec = zeros(1,N);
G_xy_vec = zeros(1,N);
nu_xy_vec = zeros(1,N);
for i=1:N
S_glob_value = S_glob(theta_vec(i));
E_x_vec(i) = 1/S_glob_value(1,1);
E_y_vec(i) = 1/S_glob_value(2,2);
G_xy_vec(i) = 1/S_glob_value(3,3);
nu_xy_vec(i) = - S_glob_value(2,1)/S_glob_value(1,1);

end

figure(1)
plot(convang(theta_vec,'rad','deg'),E_x_vec/E_2,'LineStyle',':','Color',"b","LineWidth",2,"Marker","o");
hold on;
grid on
plot(convang(theta_vec,'rad','deg'),E_y_vec/E_2,'LineStyle','--','Color',"r","LineWidth",2,"Marker","o");
title('Moduli elastici della lamina');
xlabel('$\theta (deg)$','Interpreter','latex','FontSize',12);
ylabel('GPa')
lgd = legend('$E_{x} / E_{2}$','$E_{y} / E_{2}$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;



figure(2)
plot(convang(theta_vec,'rad','deg'),G_xy_vec/E_2,'LineStyle','-.','Color',"b","LineWidth",2,"Marker","o");
hold on;
grid on
plot(convang(theta_vec,'rad','deg'),nu_xy_vec,'LineStyle','-','Color',"y","LineWidth",2,"Marker","o");
title('Moduli elastici della lamina');
xlabel('$\theta (deg)$','Interpreter','latex','FontSize',12);
lgd = legend('$G_{x,y} / E_{2}$','$\nu_{x,y}$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;

