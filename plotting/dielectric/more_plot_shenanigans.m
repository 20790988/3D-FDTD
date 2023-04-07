clear
N_x = 200;
N_y = N_x;
N_z = N_x;
delta_x = 1.5e-3;

index = delta_x*(-N_y/2+1:N_y/2)+0.5*delta_x;

load("frees_N100.mat");
H_sim_50 = E_tot_line;
load("frees_N180.mat");
H_sim_100 = E_tot_line;


H_cst_50 = readmatrix("t1.txt");
H_cst_100 = readmatrix("t2.txt");


figure(1);
clf
hold on
plot(index,H_sim_50,'r');
plot(H_cst_50(:,1)/1e3,H_cst_50(:,2),'b')
xline((120-100)*delta_x);
xline((160-100)*delta_x);
ylim([0 0.08])

xlim([-0.15, 0.15])
xlabel('x (m)')
ylabel('|E_{tot}| (V/m)');
legend('Own','CST')
hold off
grid on


figure(2);
clf
hold on
plot(index,H_sim_100,'r');
plot(H_cst_100(:,1)/1e3,H_cst_100(:,2),'b')
xline((120-100)*delta_x);
xline((160-100)*delta_x);
ylim([0 0.08])

xlim([-0.15, 0.15])
xlabel('x (m)')
ylabel('|E_{tot}| (V/m)');
legend('Own','CST')
hold off
grid on