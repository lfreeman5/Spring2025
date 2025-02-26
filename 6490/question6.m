load('question2.mat')
close all;

%%% Part A %%%
% Spatial FFTs
FS_x1 = 1/x1(2);
FS_x2 = 1/x2(2);
E11_k = fft(R11_x1);
E11_k = abs(E11_k/size(x1,2)).^2;
E11_k = E11_k(1:ceil(size(x1,2)/2));
f_E11_k = FS_x1*(0:size(x1,2)/2)/size(x1,2);
f_E11_k_norm = f_E11_k * eta;
E11_k_norm = E11_k / (nu^5 * diss)^0.25;

E22_k = fft(R22_x2);
E22_k = abs(E22_k/size(x2,1)).^2;
E22_k = E22_k(1:ceil(size(x2,1)/2));
f_E22_k = FS_x2*(0:size(x2,1)/2)/size(x2,1);
E22_k_norm = E22_k / (nu^5 * diss)^0.25;
f_E22_k_norm = f_E22_k * eta;

% Temporal FFTs
t_sampled = t(t_indices);
t_sampled = squeeze(t_sampled(1,1,1,:));
% dt_sampled = mean(diff(t_sampled));
% dt_sampled = t(2);
FS_t = 1/dt_sampled;
E11_w = fft(R11_t);
E11_w = abs(E11_w/length(t_sampled)).^2; 
E11_w = E11_w(1:ceil(length(t_sampled)/2)); 
f_E11_w = FS_t * (0:ceil(length(t_sampled)/2)-1) / length(t_sampled);
f_E11_w_norm = f_E11_w * eta / U_T;
E11_w_norm = E11_w / (nu^5 * diss)^0.25;

figure;
hold on
plot(f_E11_k_norm, E11_k_norm, 'k-', 'LineWidth', 2);
plot(f_E22_k_norm, E22_k_norm, 'b-', 'LineWidth', 2);
plot(f_E11_w_norm, E11_w_norm, 'g-', 'LineWidth', 2);
hold off
xlabel('Normalized wavenumber ($\kappa \eta$ and $\omega \eta / U_T$)', 'Interpreter', 'latex');
ylabel('Normalized Energy Content ($E / (\nu^5 \epsilon)^{1/4}$)', 'Interpreter', 'latex');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');