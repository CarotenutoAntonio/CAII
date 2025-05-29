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
    Q_12 Q_22 0;  %
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


%% DIVERSE SEQUENZE di laminazione 
%% valuto effetto variazione di % di +- 45° e 90°
%% consideriamo solo laminato simmetrico (non tengo conto particolare posizionamento)
%% fisso numero lamine a 20
N_lam = 20;
seq_0 = zeros(1,N_lam);
seq_pm45 = [45, -45];
seq_pm45 = [seq_pm45, seq_pm45, seq_pm45, seq_pm45, seq_pm45, seq_pm45 ...
            seq_pm45, seq_pm45, seq_pm45, seq_pm45];
seq_pm45 = convang(seq_pm45,'deg','rad');
seq_90 = convang(90,'deg','rad')*ones(1,N_lam);
N_perc = N_lam/2;
perc_vec=linspace(0,100,N_perc+1);
E_x_mat = nan(N_perc);
G_xy_mat = nan(N_perc);
nu_xy_mat = nan(N_perc);
E_x_1vec = nan(N_perc,1);
G_xy_1vec = nan(N_perc,1);
nu_xy_1vec = nan(N_perc,1);
E_x_2vec = nan(N_perc,1);
G_xy_2vec = nan(N_perc,1);
nu_xy_2vec = nan(N_perc,1);


for i=1:(N_perc-1)
    for j=1:(N_perc-i)
        %% al variare di i varia la percentuale di lamine a 0°
        %% al variare di j varia la percentuale di +-45°
        if (i+j)<N_perc
            seq_theta=[seq_0(1:(2*i)),seq_pm45(1:2*j), seq_90(1:(N_lam-2*(i+j)))];
        elseif (i<N_perc) && (i+j==N_perc) 
                seq_theta=[seq_0(1:(2*i)),seq_pm45(1:2*j)];
        end
        seq_theta_rad = [seq_theta, fliplr(seq_theta)];

        
        N = length(seq_theta_rad);
        z_vec = t * linspace((-N/2),(N/2),N+1);
        %% costruzione matrici A, B, D
        A = zeros(3);
        B = zeros(3);
        D = zeros(3);   
        for k=1:N
            A = A + Q_glob(seq_theta_rad(k)) * (z_vec(k+1) - z_vec(k));
            B = B + Q_glob(seq_theta_rad(k)) * ((z_vec(k+1))^2 - (z_vec(k))^2)/2;
            D = D + Q_glob(seq_theta_rad(k)) * ((z_vec(k+1))^3 - (z_vec(k))^3)/3;
        end
        ABBD = [A B; B D];
        ABBD_inv = inv(ABBD);  
        %% moduli elastici effettivi      
        E_x = 1/(ABBD_inv(1,1)*t*N*6.89);
        G_xy = 1/(ABBD_inv(3,3)*t*N);
        ni_xy = - ABBD_inv(1,2) / ABBD_inv(1,1);
        E_x_mat(i,j) = E_x;
        G_xy_mat(i,j) = G_xy;
        nu_xy_mat(i,j) = ni_xy;
        %% moduli di bending
        E_bx = 12/(ABBD_inv(4,4)*(t*N)^3);
        E_by = 12/(ABBD_inv(5,5)*(t*N)^3);
    end
end

%% caso non ci sono lamine a 0°
 for j=1:N_perc
        %% al variare di i varia la percentuale di lamine a 0°
        %% al variare di j varia la percentuale di +-45°
        if (j<N_perc)

            seq_theta=[seq_pm45(1:2*j), seq_90(1:(N_lam-2*(j)))];

        elseif (j==N_perc) 
                seq_theta=seq_pm45(1:2*j);
        end

        seq_theta_rad = [seq_theta, fliplr(seq_theta)];

        
        N = length(seq_theta_rad);
        z_vec = t * linspace((-N/2),(N/2),N+1);
        %% costruzione matrici A, B, D
        A = zeros(3);
        B = zeros(3);
        D = zeros(3);
        
        for k=1:N
            A = A + Q_glob(seq_theta_rad(k)) * (z_vec(k+1) - z_vec(k));
            B = B + Q_glob(seq_theta_rad(k)) * ((z_vec(k+1))^2 - (z_vec(k))^2)/2;
            D = D + Q_glob(seq_theta_rad(k)) * ((z_vec(k+1))^3 - (z_vec(k))^3)/3;
        end
        ABBD = [A B; B D];
        ABBD_inv = inv(ABBD);
        
        %% moduli elastici effettivi
        
        E_x = 1/(ABBD_inv(1,1)*t*N*6.89);
        G_xy = 1/(ABBD_inv(3,3)*t*N);
        ni_xy = - ABBD_inv(1,2) / ABBD_inv(1,1);

        E_x_1vec(j) = E_x;

        G_xy_1vec(j) = G_xy;

        nu_xy_1vec(j) = ni_xy;


        %% moduli di bending
        
        E_bx = 12/(ABBD_inv(4,4)*(t*N)^3);
        E_by = 12/(ABBD_inv(5,5)*(t*N)^3);
 end


 %% caso non ci sono lamine a +-45°
 for j=1:N_perc

        if (j<N_perc)

            seq_theta=[seq_0(1:2*j), seq_90(1:(N_lam-2*(j)))];

        elseif (j==N_perc) 
                seq_theta=seq_0(1:2*j);
        end

        seq_theta_rad = [seq_theta, fliplr(seq_theta)];

        
        N = length(seq_theta_rad);
        z_vec = t * linspace((-N/2),(N/2),N+1);
        %% costruzione matrici A, B, D
        A = zeros(3);
        B = zeros(3);
        D = zeros(3);
        
        for k=1:N
            A = A + Q_glob(seq_theta_rad(k)) * (z_vec(k+1) - z_vec(k));
            B = B + Q_glob(seq_theta_rad(k)) * ((z_vec(k+1))^2 - (z_vec(k))^2)/2;
            D = D + Q_glob(seq_theta_rad(k)) * ((z_vec(k+1))^3 - (z_vec(k))^3)/3;
        end
        ABBD = [A B; B D];
        ABBD_inv = inv(ABBD);
        
        %% moduli elastici effettivi
        
        E_x = 1/(ABBD_inv(1,1)*t*N*6.89);
        G_xy = 1/(ABBD_inv(3,3)*t*N);
        ni_xy = - ABBD_inv(1,2) / ABBD_inv(1,1);

        E_x_2vec(j) = E_x;

        G_xy_2vec(j) = G_xy;

        nu_xy_2vec(j) = ni_xy;


        %% moduli di bending
        
        E_bx = 12/(ABBD_inv(4,4)*(t*N)^3);
        E_by = 12/(ABBD_inv(5,5)*(t*N)^3);
 end

%% caso tutte lamine a 90°
        seq_theta_rad = seq_90;
        N = length(seq_theta_rad);
        z_vec = t * linspace((-N/2),(N/2),N+1);
        z_mean= z_vec(2:end)-t/2;
        %% costruzione matrici A, B, D
        A = zeros(3);
        B = zeros(3);
        D = zeros(3);
        
        for k=1:N
            A = A + Q_glob(seq_theta_rad(k)) * (z_vec(k+1) - z_vec(k));
            B = B + Q_glob(seq_theta_rad(k)) * ((z_vec(k+1))^2 - (z_vec(k))^2)/2;
            D = D + Q_glob(seq_theta_rad(k)) * ((z_vec(k+1))^3 - (z_vec(k))^3)/3;
        end
        ABBD = [A B; B D];
        ABBD_inv = inv(ABBD);
        
        %% moduli elastici effettivi
        
        E_x = 1/(ABBD_inv(1,1)*t*N*6.89);
        E_y = 1/(ABBD_inv(2,2)*t*N);
        G_xy = 1/(ABBD_inv(3,3)*t*N);
        ni_xy = - ABBD_inv(1,2) / ABBD_inv(1,1);





 E_x_2vec = [E_x, E_x_2vec'];
 
 E_x_mat = [E_x_1vec'; E_x_mat];

 E_x_mat = [E_x_2vec(:), E_x_mat];


  G_xy_2vec = [G_xy, G_xy_2vec'];
 
 G_xy_mat = [G_xy_1vec'; G_xy_mat];

 G_xy_mat = [G_xy_2vec(:), G_xy_mat];


  nu_xy_2vec = [ni_xy, nu_xy_2vec'];
 
 nu_xy_mat = [nu_xy_1vec'; nu_xy_mat];

 nu_xy_mat = [nu_xy_2vec(:), nu_xy_mat];



 %% grafica carpet plot
   figure(1)
   v=linspace(1,N_perc+1,N_perc+1);
   v_s=fliplr(v);

    plot(perc_vec,E_x_mat(1,:),"Color",[0 0.4470 0.7410]);
    hold on,
    grid on,
    plot(perc_vec,E_x_mat(3,:),"Color",[0.8500 0.3250 0.0980]);
    plot(perc_vec,E_x_mat(5,:),"Color",[0.9290 0.6940 0.1250]);
    plot(perc_vec,E_x_mat(7,:),"Color",[0.4940 0.1840 0.5560]);
    plot(perc_vec,E_x_mat(9,:),"Color",[0.4660 0.6740 0.1880]);
    plot(perc_vec,E_x_mat(11,:),"Color",[0.3010 0.7450 0.9330]);
    %plot([perc_vec(3)],[E_x_mat(3,3),E_x_mat(5,5),E_x_mat(7,7),E_x_mat(9,9),E_x_mat(11,11)])
    plot(flip(perc_vec),diag(fliplr(E_x_mat)));
    title('CARPET PLOT');
ylabel('$E_{x}(MSI)$','Interpreter','latex','FontSize',12);
xlabel(' percentuale $\pm45^\circ$ ','Interpreter','latex','FontSize',12);
lgd = legend(' 0% 0° ',' 20% 0° ',' 40% 0° ',' 60% 0° ','80% 0°','100% 0°');
%lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;

figure(2)
    plot(perc_vec,nu_xy_mat(1,:),"Color",[0 0.4470 0.7410]);
    hold on,
    grid on,
    plot(perc_vec,nu_xy_mat(3,:),"Color",[0.8500 0.3250 0.0980]);
    plot(perc_vec,nu_xy_mat(5,:),"Color",[0.9290 0.6940 0.1250]);
    plot(perc_vec,nu_xy_mat(7,:),"Color",[0.4940 0.1840 0.5560]);
    plot(perc_vec,nu_xy_mat(9,:),"Color",[0.4660 0.6740 0.1880]);
    plot(perc_vec,nu_xy_mat(11,:),"Color",[0.3010 0.7450 0.9330]);
    plot(flip(perc_vec),diag(fliplr(nu_xy_mat)));
    title('CARPET PLOT');
ylabel('$nu_{xy}(MSI)$','Interpreter','latex','FontSize',12);
xlabel(' percentuale $\pm45^\circ$ ','Interpreter','latex','FontSize',12);
lgd = legend(' 0% 0\circ ',' 20% 0° ',' 40% 0° ',' 60% 0° ','80% 0°','100% 0°');
%lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;


figure(3)
    plot(perc_vec,G_xy_mat(1,:),"Color",[0 0.4470 0.7410]);

    title('CARPET PLOT');
ylabel('$G_{xy}(MSI)$','Interpreter','latex','FontSize',12);
xlabel(' percentuale $\pm45^\circ$ ','Interpreter','latex','FontSize',12);
