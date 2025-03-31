clc; close all;

%% Part A %%
% Calculate fluctuating velocity components
t_rand = 180;
u1_prime = u1 - mean(u1, 4);
u2_prime = u2 - mean(u2, 4);
w3 = zeros(size(u1_prime(:,:,1,t_rand)));
w3(2:end-1,2:end-1)=((u2(3:end,2:end-1,1,t_rand)-u2(1:end-2,2:end-1,1,t_rand))/(2*(x1(2)-x1(1)))) - ...
    (u1(2:end-1,3:end,1,t_rand)-u1(2:end-1,1:end-2,1,t_rand))/(2*(x2(2)-x2(1)));


figure
imagesc(x1, x2, u1_prime(:,:,1,t_rand))
set(gca, 'ydir', 'normal')
xlabel('x(m)')
ylabel('y(m)')
colormap(jet)
cb = colorbar;
ylabel(cb,'U1'' (m/s)')

figure
imagesc(x1, x2, u2_prime(:,:,1,t_rand))
set(gca, 'ydir', 'normal')
xlabel('x(m)')
ylabel('y(m)')
colormap(jet)
cb = colorbar;
ylabel(cb,'U2'' (m/s)')

figure
imagesc(x1, x2, w3)
set(gca, 'ydir', 'normal')
xlabel('x(m)')
ylabel('y(m)')
colormap(jet)
cb = colorbar;
ylabel(cb,'wz (1/s)')

%% Part B %%
u1_bar = mean(u1,4);
u2_bar = mean(u2,4);
ubar_mag = sqrt(u1_bar.^2 + 2 * u2_bar.^2);
u1_rms = sqrt(mean(u1.^2, 4));
u2_rms = sqrt(mean(u2.^2, 4));
u_rms = (u1_rms+2*u2_rms)/3;

figure
quiver(x1,x2,u1_bar,u2_bar,2)
xlabel('x (m)')
ylabel('y (m)')
axis tight

figure
imagesc(x1, x2, ubar_mag)
set(gca, 'ydir', 'normal')
xlabel('x (m)')
ylabel('y (m)')
colormap(jet)
cb = colorbar;
ylabel(cb,'mag(u bar) (m/s)')

figure
imagesc(x1, x2, u1_rms)
set(gca, 'ydir', 'normal')
xlabel('x (m)')
ylabel('y (m)')
colormap(jet)
cb = colorbar;
ylabel(cb,'rms(u1) (m/s)')

figure
imagesc(x1, x2, u2_rms)
set(gca, 'ydir', 'normal')
xlabel('x (m)')
ylabel('y (m)')
colormap(jet)
cb = colorbar;
ylabel(cb,'rms(u2) (m/s)')

figure
imagesc(x1, x2, u_rms)
set(gca, 'ydir', 'normal')
xlabel('x (m)')
ylabel('y (m)')
colormap(jet)
cb = colorbar;
ylabel(cb,'rms(u) (m/s)')

%% Part C %% 
u1_rms_avg = mean(u1_rms, 'all');
u2_rms_avg = mean(u2_rms, 'all');
u_rms_avg = mean(u_rms, 'all');
u_rms_ratio = u1_rms_avg/u2_rms_avg;
k = 0.5 * (mean(u1_prime.^2,4) + 2*mean(u2_prime.^2,4));
k_avg = mean(k, 'all');

%% Part D %%
mff = ubar_mag./u_rms;
mff_avg = mean(mff,'all');

figure
imagesc(x1, x2, mff)
set(gca, 'ydir', 'normal')
xlabel('x (m)')
ylabel('y (m)')
colormap(jet)
cb = colorbar;
ylabel(cb,'Mean Flow Factor (dimensionless)')

%% Part E %%
k_spatial = 0.5 * (mean(u1_prime.^2,[1,2]) + 2*mean(u2_prime.^2,[1,2]));

figure
plot(squeeze(t), squeeze(k_spatial));
xlabel('Time (s)')
ylabel('Spatially Averaged Kinetic Energy (J/kg)')
