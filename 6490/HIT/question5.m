close all;

%%% Part A %%%
eta = (nu^3 / diss)^0.25;
tau = sqrt(nu/diss);
u_eta = (nu*diss)^0.25;
a_eta = (diss^3 /nu)^0.25;

%%% Part B %%%
taylor_a = sqrt(15*nu*u_rms_avg^2/diss);
N = 3;
k1 = sum((1 - rho11_x1(1:N)) .* x1(1:N).^2) / sum(x1(1:N).^4);
k2 = sum((1 - rho22_x2(1:N)) .* x2(1:N)'.^2) / sum(x2(1:N).^4);
taylor_b1 = sqrt(1/k1);
taylor_b2 = sqrt(1/k2);
