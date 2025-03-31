% clearvars -except x1 x2 t u1 u2 u1_prime u2_prime
close all
clc

if exist('question3.mat', 'file') == 2
    load('question3.mat')
else
    %%% Part A %%%
    D11_x1 = zeros(1,size(x1,2)-1);
    D22_x1 = zeros(1,size(x1,2)-1);
    for offset = 1:size(x1,2)-1
        D11_x1(offset) = mean((u1_prime(:,1:end-offset,1,:)- ...
                         u1_prime(:,1+offset:end,1,:)).^2,'all');
        D22_x1(offset) = mean((u2_prime(:,1:end-offset,1,:)- ...
                         u2_prime(:,1+offset:end,1,:)).^2,'all');
    end
    
    D11_x2 = zeros(1,size(x2,1)-1);
    D22_x2 = zeros(1,size(x2,1)-1);
    for offset = 1:size(x2,1)-1
        D11_x2(offset) = mean((u1_prime(1:end-offset,:,1,:)- ...
                         u1_prime(offset+1:end,:,1,:)).^2, 'all');
        D22_x2(offset) = mean((u2_prime(1:end-offset,:,1,:)- ...
                         u2_prime(offset+1:end,:,1,:)).^2, 'all');
    end

    save('question3.mat', 'D11_x1', 'D11_x2', 'D22_x1', 'D22_x2');
end

%%% Plot part A %%%
figure;
hold on;
plot(x1(2:end), D11_x1, 'r-', 'LineWidth', 2); % Red solid line
plot(x2(2:end), D22_x2, 'b--', 'LineWidth', 2); % Blue dashed line
func = @(x) x.^2 * D11_x1(1)/(x1(2)^2); 
plot(x1(2:end), func(x1(2:end)), 'g--', 'LineWidth', 2); % Red solid line
func = @(x) x.^(2/3) * D11_x1(20)/(x1(21)^(2/3));
plot(x1(2:end), func(x1(2:end)), 'c--', 'LineWidth', 2); % Red solid line


plot(x2(2:end), D11_x2, 'r--', 'LineWidth', 2); % Red dashed line
plot(x1(2:end), D22_x1, 'b-', 'LineWidth', 2); % Blue solid line
func = @(x) x.^2 * D11_x2(1)/(x2(2)^2); 
plot(x2(2:end), func(x2(2:end)), 'g-', 'LineWidth', 2); % Red solid line
func = @(x) x.^(2/3) * D11_x2(20)/(x2(21)^(2/3));
plot(x2(2:end), func(x2(2:end)), 'c-', 'LineWidth', 2); % Red solid line

hold off;
xlabel('l (m)');
ylabel('D');
legend({'D_{11}(x_1)', 'D_{22}(x_1)', 'x^{2} scaling for D_{LL}', 'x^{2/3} scaling for D_{LL}', 'D_{11}(x_2)', 'D_{22}(x_2)', 'x^{2} scaling for D_{LL}', 'x^{2/3} scaling for D_{LL}'}, 'Location', 'Best');
title('2nd order structure functions');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');