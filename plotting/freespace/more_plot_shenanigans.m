clear
N_x = 100;
N_y = N_x;
N_z = N_x;
delta_x = 3e-3;

index = delta_x*(-N_y/2+1:N_y/2)+0.5*delta_x;

load("frees_N50.mat");
H_sim_50 = H_tot_line;
load("frees_N100.mat");
H_sim_100 = H_tot_line;
load("frees_N200.mat");
H_sim_150 = H_tot_line;

H_cst_50 = readmatrix("t1.txt");
H_cst_100 = readmatrix("t2.txt");
H_cst_150 = readmatrix("t3.txt");


figure(1);
clf
hold on
plot(index,H_sim_50,'r');
plot(H_cst_50(:,1)/1e3,H_cst_50(:,2),'b')
ylim([0 8e-5])

xlim([-0.15, 0.15])
xlabel('x (m)')
ylabel('|H_{tot}| (A/m)');
legend('Own','CST')
hold off
grid on

figure(3);
clf
hold on
plot(index,H_sim_150,'r');
plot(H_cst_150(:,1)/1e3,H_cst_150(:,2),'b')
ylim([0 8e-5])

xlim([-0.15, 0.15])
xlabel('x (m)')
ylabel('|H_{tot}| (A/m)');
legend('Own','CST')
hold off
grid on