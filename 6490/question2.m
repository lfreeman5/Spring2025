% clearvars -except x1 x2 t u1 u2 u1_prime u2_prime
close all
clc

if exist('question2.mat', 'file') == 2
    load('question2.mat')
else
    %%% Part A %%%
    % R11_x1 = zeros(1,size(x1,2)-1);
    % R22_x1 = zeros(1,size(x1,2)-1);
    % for offset = 1:size(x1,2)-1
    %     R11_x1(offset) = mean(u1_prime(:,1:end-offset,1,:).* ...
    %                      u1_prime(:,1+offset:end,1,:),'all');
    %     R22_x1(offset) = mean(u2_prime(:,1:end-offset,1,:).* ...
    %                      u2_prime(:,1+offset:end,1,:),'all');
    % end
    % 
    % R11_x2 = zeros(1,size(x2,1)-1);
    % R22_x2 = zeros(1,size(x2,1)-1);
    % for offset = 1:size(x2,1)-1
    %     R11_x2(offset) = mean(u1_prime(1:end-offset,:,1,:).* ...
    %                      u1_prime(offset+1:end,:,1,:), 'all');
    %     R22_x2(offset) = mean(u2_prime(1:end-offset,:,1,:).* ...
    %                      u2_prime(offset+1:end,:,1,:), 'all');
    % end
    % 
    % rho11_x1 = R11_x1/mean(u1_prime.^2,'all');
    % rho22_x1 = R22_x1/mean(u2_prime.^2,'all');
    % rho11_x2 = R11_x2/mean(u1_prime.^2,'all');
    % rho22_x2 = R22_x2/mean(u2_prime.^2,'all');
    
    %%% Part B %%%
    % Nothing here for now
    
    %%% Part C %%%
    t_indices = 1:5:15250;
    % t_indices = [t_indices, 1001:10:3000];
    % t_indices = [t_indices, 3001:100:size(t,4)-1];
    R11_t = zeros(1, length(t_indices));
    R22_t = zeros(1, length(t_indices));
    for i = 1:length(t_indices)
        offset = t_indices(i);
        disp(offset)
        R11_t(i) = mean(u1_prime(:,:,1,1:end-offset).* ...
                         u1_prime(:,:,1,offset+1:end), 'all');
        R22_t(i) = mean(u2_prime(:,:,1,1:end-offset).* ...
                         u2_prime(:,:,1,offset+1:end), 'all');
    end
    rho11_t = R11_t/mean(u1_prime.^2,'all');
    rho22_t = R22_t/mean(u2_prime.^2,'all');
    
end

%%% Plot part A %%%
% figure;
% hold on;
% plot(x1(2:end), rho11_x1, 'r-', 'LineWidth', 2); % Red solid line
% plot(x1(2:end), rho22_x1, 'b-', 'LineWidth', 2); % Blue solid line
% plot(x2(2:end), rho11_x2, 'r--', 'LineWidth', 2); % Red dashed line
% plot(x2(2:end), rho22_x2, 'b--', 'LineWidth', 2); % Blue dashed line
% hold off;
% xlabel('l (m)');
% ylabel('Correlation');
% legend({'\rho_{11}(x_1)', '\rho_{22}(x_1)', '\rho_{11}(x_2)', '\rho_{22}(x_2)'}, 'Location', 'Best');
% title('Velocity Correlations');
% grid on;
% 
% %%% Calculate part B %%%
% trim_start = 11; trim_end_x1 = length(x1)-5; trim_end_x2 = length(x2)-5;
% x1_t = x1(trim_start:trim_end_x1)'; % Transpose to column
% rho11_x1_t = rho11_x1(trim_start:trim_end_x1)'; 
% rho22_x1_t = rho22_x1(trim_start:trim_end_x1)'; 
% x2_t = x2(trim_start:trim_end_x2); % Already column, no transpose needed
% rho11_x2_t = rho11_x2(trim_start:trim_end_x2)'; 
% rho22_x2_t = rho22_x2(trim_start:trim_end_x2)';
% f_exp = fittype('A*exp(B*x)', 'independent', 'x', 'coefficients', {'A', 'B'});
% fit_rho11_x1 = fit(x1_t, rho11_x1_t, f_exp); 
% fit_rho22_x1 = fit(x1_t, rho22_x1_t, f_exp);
% fit_rho11_x2 = fit(x2_t, rho11_x2_t, f_exp); 
% fit_rho22_x2 = fit(x2_t, rho22_x2_t, f_exp);
% 
% rho_target = 0.4;
% L_11_1 = log(rho_target / fit_rho11_x1.A) / fit_rho11_x1.B;
% L_22_1 = log(rho_target / fit_rho22_x1.A) / fit_rho22_x1.B;
% L_11_2 = log(rho_target / fit_rho11_x2.A) / fit_rho11_x2.B;
% L_22_2 = log(rho_target / fit_rho22_x2.A) / fit_rho22_x2.B;
% L_LL = (L_11_1+2*L_22_2)/3;
% 
% 
% %%% Plot part C %%%
% figure;
% plot(squeeze(t(1,1,1,t_indices)), rho11_t, 'r-', 'LineWidth', 2); % R11_t in red
% hold on;
% plot(squeeze(t(1,1,1,t_indices)), rho22_t, 'b-', 'LineWidth', 2); % R22_t in blue
% xlabel('t (s)');
% ylabel('Temporal Correlation');
% title('Temporal Correlation \rho11_t and \rho22_t vs Time');
% legend({'\rho11_t', '\rho22_t'}, 'Location', 'Best');
% grid on;
% hold off;
% 
% figure;
% filtered_indices = t_indices(t_indices < 10000);
% filtered_r11_t = rho11_t(t_indices<10000);
% filtered_r22_t = rho22_t(t_indices<10000);
% plot(squeeze(t(1,1,1,filtered_indices)), filtered_r11_t, 'r-', 'LineWidth', 2);
% hold on;
% plot(squeeze(t(1,1,1,filtered_indices)), filtered_r22_t, 'b-', 'LineWidth', 2);
% xlabel('t (s)');
% ylabel('Temporal Correlation');
% title('Temporal Correlation \rho11_t and \rho22_t vs Time');
% legend({'\rho11_t', '\rho22_t'}, 'Location', 'Best');
% grid on;
% hold off;
% 
% %%% Calculate Part D %%%
% T_11 = t(t_indices(find(rho11_t < 0.4, 1)));
% T_22 = t(t_indices(find(rho22_t < 0.4, 1)));
% T_LL = (T_11+2*T_22)/3.;
% 
% %%% Calculate Part E %%%
% U_T = L_LL/T_LL;