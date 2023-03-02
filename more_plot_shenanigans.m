clear
N_x = 100;
N_y = N_x;
N_z = N_x;
delta_x = 3e-3;

index = delta_x*(-N_y/2+1:N_y/2)+0.5*delta_x;

load("step50.mat");
H_sim_50 = H_tot(N_x/2,1:N_y,N_z/2);
load("step100.mat");
H_sim_100 = H_tot(N_x/2,1:N_y,N_z/2);
load("step150.mat");
H_sim_150 = H_tot(N_x/2,1:N_y,N_z/2);

H_cst_50 = readmatrix("t1.txt");
H_cst_100 = readmatrix("t2.txt");
H_cst_150 = readmatrix("t3.txt");

figure(1);
clf
hold on
plot(index,H_sim_50);
plot(H_cst_50(:,1)/1e3,H_cst_50(:,2))
ylim([0 2e-4])
hold off
grid on

figure(2);
clf
hold on
plot(index,H_sim_100);
plot(H_cst_100(:,1)/1e3,H_cst_100(:,2))
ylim([0 2e-4])
hold off
grid on

figure(3);
clf
hold on
plot(index,H_sim_150);
plot(H_cst_150(:,1)/1e3,H_cst_150(:,2))
ylim([0 2e-4])
hold off
grid on