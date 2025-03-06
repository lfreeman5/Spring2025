close all;
clc;

nu = 1.46*10^(-5); %m^2/s
dx = x1(2);
dy = x2(2);

%%% Part A %%%
du1p_dx1 = (u1_prime(2:end-1,3:end,:,:) - u1_prime(2:end-1,1:end-2,:,:))/(2*dx);
du1p_dx2 = (u1_prime(3:end,2:end-1,:,:) - u1_prime(1:end-2,2:end-1,:,:))/(2*dy);
du2p_dx1 = (u2_prime(2:end-1,3:end,:,:) - u2_prime(2:end-1,1:end-2,:,:))/(2*dx);
du2p_dx2 = (u2_prime(3:end,2:end-1,:,:) - u2_prime(1:end-2,2:end-1,:,:))/(2*dy);
pseudo_diss_a = nu * (-(du1p_dx1.^2) + 2*(du1p_dx2.^2) + 2*(du2p_dx1.^2)+8*(du2p_dx2.^2));
diss_a_spatial = mean(pseudo_diss_a,4);
diss_a = mean(diss_a_spatial, 'all');

%%% Calculate Part B %%%
C2 = 2.0;
D11_x1_compensated = (D11_x1/C2).^1.5 ./ x1(2:end);
D22_x2_compensated = (D22_x2/C2).^1.5 ./ x2(2:end)';
D11_x2_compensated = (D11_x2*0.75/C2).^1.5 ./ x2(2:end)';
D22_x1_compensated = (D22_x1*0.75/C2).^1.5 ./ x1(2:end);

x1_indices = find(x1(1:end-1) >= 0.012 & x1(1:end-1) <= 0.056); 
x2_indices = find(x2(1:end-1) >= 0.012 & x2(1:end-1) <= 0.056); 

diss_b = mean([mean(D11_x1_compensated(x1_indices)), ...
                      mean(D22_x2_compensated(x2_indices)), ...
                      mean(D11_x2_compensated(x2_indices)), ...
                      mean(D22_x1_compensated(x1_indices))]);

%%% Calculate Part C %%%
C_eps = 0.5;
diss_c = C_eps*u_rms_avg^3 / L_LL;

%%% Calculate Part D %%%
du1b_dx1 = (u1_bar(2:end-1,3:end,:,:) - u1_bar(2:end-1,1:end-2,:,:))/(2*dx);
du1b_dx2 = (u1_bar(3:end,2:end-1,:,:) - u1_bar(1:end-2,2:end-1,:,:))/(2*dy);
du2b_dx1 = (u2_bar(2:end-1,3:end,:,:) - u2_bar(2:end-1,1:end-2,:,:))/(2*dx);
du2b_dx2 = (u2_bar(3:end,2:end-1,:,:) - u2_bar(1:end-2,2:end-1,:,:))/(2*dy);
u1_prime_trimmed = u1_prime(2:end-1, 2:end-1, :, :);
u2_prime_trimmed = u2_prime(2:end-1, 2:end-1, :, :);
P = -mean(u1_prime_trimmed.^2,4).*du1b_dx1 ...
    - 2*mean(u1_prime_trimmed.*u2_prime_trimmed,4).*(du1b_dx2 + du2b_dx1) ...
    - 2*mean(u2_prime_trimmed.^2,4).*du2b_dx2;
diss_d = mean(P,'all');

%%% Plotting %%%
figure;
imagesc(x1, x2, pseudo_diss_a(:,:,t_rand)'); 
set(gca, 'YDir', 'normal'); 
colormap jet; colorbar;

xlabel('$x_1$ (m)', 'Interpreter', 'latex');
ylabel('$x_2$ (m)', 'Interpreter', 'latex');
ylabel(colorbar, '$\tilde{\epsilon} $(m$^2$/s$^3$)', 'Interpreter', 'latex');
title('Scalar Color Plot of $\tilde{\epsilon}$ at Timestep 180', 'Interpreter', 'latex');

figure;
imagesc(x1(3:end-2), x2(3:end-2), diss_a_spatial(3:end-2, 3:end-2)'); 
set(gca, 'YDir', 'normal'); 
colormap jet; colorbar;
clim([0, 1.1*max(diss_a_spatial(:))]);
xlabel('$x_1$ (m)', 'Interpreter', 'latex');
ylabel('$x_2$ (m)', 'Interpreter', 'latex');
ylabel(colorbar, '$\epsilon $(m$^2$/s$^3$)', 'Interpreter', 'latex');
title('Scalar Color Plot of $\epsilon (x,y)$', 'Interpreter', 'latex');

% Should compensated structure functions be log-scaled?
figure;
hold on;
plot(x1(2:end), D11_x1_compensated, 'r-', 'LineWidth', 2); % Red solid line
plot(x2(2:end), D22_x2_compensated, 'b--', 'LineWidth', 2); % Blue dashed line
plot(x2(2:end), D11_x2_compensated, 'r--', 'LineWidth', 2); % Red dashed line
plot(x1(2:end), D22_x1_compensated, 'b-', 'LineWidth', 2); % Blue solid line

hold off;
xlabel('l (m)');
ylabel('Compensated 2nd-order structure function, m^{2}/s^{2}');
legend({'D_{11}(x_1) [Compensated]', 'D_{22}(x_1) [Compensated]', 'D_{11}(x_2) [Compensated]', 'D_{22}(x_2) [Compensated]'}, 'Location', 'Best');
title('Compensated 2nd order structure functions');
grid on;
% set(gca, 'XScale', 'log', 'YScale', 'log');

% Choose dissipation
diss = diss_b;