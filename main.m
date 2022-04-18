%% clean environment
clear
clc
close all
cvx_clear;


%% setup system

sampling_period              = 60*5;                                            % sec
orbital_radius               = 35622 + 6378.1;                                  % km
gravitational_constant       = 6.673e-20;                                       % km^3 kg^-1 sec^-2
celestial_mass               = 5.9472e24;                                       % kg
gravitational_body           = gravitational_constant * celestial_mass;         % km^3 sec^-2
orbit_ang_vel                = sqrt(gravitational_body / orbital_radius^3);     % rad sec^-2

% Continuous-time LTI CWH unforced dynamics e^{At}
e_power_At = @(t) [ 
    4 - 3 * cos(orbit_ang_vel * t), 0, 0, (1/orbit_ang_vel) * sin(orbit_ang_vel * t), (2/orbit_ang_vel) * (1 - cos(orbit_ang_vel * t)), 0; 
    6 * (sin(orbit_ang_vel * t) - orbit_ang_vel * t), 1, 0, -(2/orbit_ang_vel) * (1 - cos(orbit_ang_vel * t)), (1/orbit_ang_vel) * (4*sin(orbit_ang_vel * t) - 3*orbit_ang_vel * t), 0; 
    0, 0, cos(orbit_ang_vel * t), 0, 0, (1/orbit_ang_vel) * sin(orbit_ang_vel * t); 
    3 * orbit_ang_vel * sin(orbit_ang_vel * t), 0, 0, cos(orbit_ang_vel * t), 2 * sin(orbit_ang_vel * t), 0; 
    -6 * orbit_ang_vel * (1 - cos(orbit_ang_vel * t)), 0, 0, -2 * sin(orbit_ang_vel * t), 4 * cos(orbit_ang_vel * t) - 3, 0;
    0, 0, -orbit_ang_vel * sin(orbit_ang_vel * t), 0, 0, cos(orbit_ang_vel * t);
    ];

% Discrete-time system is Phi(T_s) for sampling time T_s since the system is time-invariant
Ad = e_power_At(sampling_period);

% Impulse control
Bd = Ad*[zeros(3); eye(3)];

%% problem set up
time_horizon = 8;

% initial states
% format: x, y, z,  x., y., z.
x_0_a = [90;  -5; 0; 0; 0; 0] ; % satellite A
x_0_b = [95;   5; 0; 0; 0; 0] ; % satellite B
x_0_c = [100; -5; 0; 0; 0; 0] ; % satellite C
x_0_d = [105;  5; 0; 0; 0; 0] ; % satellite E
x_0_e = [110; -5; 0; 0; 0; 0] ; % satellite D

% target sets
% format: x, y, z, x., y., z.
target_set_a = Polyhedron('lb', [-7.5;   -10; -2.5; -0.01; -0.01; -0.01], ... 
                          'ub', [-2.5;    -5;  2.5;  0.01; 0.01;  0.01]);  
target_set_b = Polyhedron('lb', [-12.5; -2.5; -2.5; -0.01; -0.01; -0.01], ...
                          'ub', [-7.5;   2.5;  2.5;  0.01; 0.01;  0.01]);    
target_set_c = Polyhedron('lb', [-7.5;     5; -2.5; -0.01; -0.01; -0.01], ... 
                          'ub', [-2.5;    10;  2.5;  0.01; 0.01;  0.01]);   
target_set_d = Polyhedron('lb', [ 2.5;  -7.5; -2.5; -0.01; -0.01; -0.01], ...
                          'ub', [ 7.5;  -2.5;  2.5;  0.01; 0.01;  0.01]);    
target_set_e = Polyhedron('lb', [ 2.5;   2.5; -2.5; -0.01; -0.01; -0.01], ... 
                          'ub', [ 7.5;   7.5;  2.5;  0.01; 0.01;  0.01]);   
                      
n_lin_state_a = size(target_set_a.A,1);
n_lin_state_b = size(target_set_b.A,1);
n_lin_state_c = size(target_set_c.A,1);
n_lin_state_d = size(target_set_d.A,1);
n_lin_state_e = size(target_set_e.A,1);
                      
% Input space
input_space = Polyhedron('lb', [ -3; -3; -3], ... 
                          'ub', [ 3;  3;  3]);                         

input_space_A = blkdiag(input_space.A);
for i=1:(time_horizon-1)
    input_space_A = blkdiag(input_space_A, input_space.A);
end

input_space_b = repmat(input_space.b, time_horizon,1);

% collision avoid region radius
r = 5;

% safety threshold
safety_target       = 0.8; 
safety_collision    = 0.8;



%% Concat matrix maker

Ad_concat = zeros(size(Ad, 1)*time_horizon, size(Ad, 2));
Bd_concat = zeros(size(Bd, 1)*time_horizon, size(Bd, 2)*time_horizon);
Wd_concat = zeros(size(Ad, 1)*time_horizon);
for i = 0:(time_horizon-1)
    Ad_concat(size(Ad, 1)*i + [1:size(Ad, 1)], :) = Ad^(i+1);
end
for i = 0:(time_horizon-1)
    for j = 0:i
        Bd_concat(size(Bd, 1)*i + [1:size(Bd, 1)], size(Bd, 2)*j + [1:size(Bd, 2)]) = Ad^(i-j) * Bd;
        Wd_concat(size(Ad, 1)*i + [1:size(Ad, 1)], size(Ad, 2)*j + [1:size(Ad, 2)]) = Ad^(i-j);
    end
end


%% disturbance set up
mu = zeros(6,1);
cov = diag([10^-4, 10^-4, 10^-4, 5*10^-8, 5*10^-8, 5*10^-8]);
psi = cov;
for i=1:(time_horizon-1)
    psi = blkdiag(psi, cov);
end
nu = 10;

Cw_psi_Cw = Wd_concat * psi * Wd_concat';

% multiplier for target set
scaled_sigma_a_vec = diag(target_set_a.A * Cw_psi_Cw(end-5:end,end-5:end) * target_set_a.A').^(1/2);
scaled_sigma_b_vec = diag(target_set_b.A * Cw_psi_Cw(end-5:end,end-5:end) * target_set_b.A').^(1/2);
scaled_sigma_c_vec = diag(target_set_c.A * Cw_psi_Cw(end-5:end,end-5:end) * target_set_c.A').^(1/2);
scaled_sigma_d_vec = diag(target_set_d.A * Cw_psi_Cw(end-5:end,end-5:end) * target_set_d.A').^(1/2);
scaled_sigma_e_vec = diag(target_set_e.A * Cw_psi_Cw(end-5:end,end-5:end) * target_set_e.A').^(1/2);

% multiplier for collision avoidance
max_sigma = zeros(time_horizon-1,1);
for i = 1:(time_horizon-1)
    index = 4*(i-1) + (1:3); 
    max_sigma(i) = 2 * nu * max(eig(Cw_psi_Cw(index, index)));
end

% disturbance  linearization
%set up
syms t_pdf(x) beta_prime_pdf(x);

by = 5e-6; % step size
approx_to = .99995;
approx_from = .5;
tolerance = 1e-1;
linearize_from = safety_collision;

lb_approx = 1 - approx_to;
ub_approx = 1 - linearize_from;

% t dist
t_init = 0; % known quantile at p=0.5
t_pdf(x) = gamma((nu+1)/2) / (sqrt(nu*pi) * gamma(nu/2)) * (1 + x^2/nu)^(-(nu+1)/2);
[t_invcdf_m, t_invcdf_c] = quantile_approx(1, by, approx_to, approx_from, t_init, t_pdf, tolerance, linearize_from);

% beta prome dist
q = size(Ad, 1)/2;
a = 2 * q * ( q/2 + nu^2/4 - nu + q*nu/2 - 2*q + 1) / ((nu - 2)*(q + nu - 2));
b = 2 * (-q + nu^2/4 - nu/2 + q*nu/2) / (q + nu - 2);

beta_prime_med_approx = 2^(-1/a) * (log(2) - 1/3 + a) / (2^(-1/b) * (log(2) - 1/3 + b));
beta_prime_pdf(x) = x^(a-1) * (1+x)^(-a-b) / beta(a,b);
[beta_invcdf_m, beta_invcdf_c] = quantile_approx(1, by, approx_to, approx_from, beta_prime_med_approx, beta_prime_pdf, tolerance, linearize_from);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
kmax = 100;

% convergence perameters
epsilon_dc = 1e-6; % convergence in cost
epsilon_lambda = 1e-6; % convergence of sum of slack variables to zero

% cost of slack variable
tau_max = 100;
gamma = 1.2;
tau = min(tau_max * ones(kmax,1), gamma.^(0:(kmax-1))');

% storage initial cost for convergence check
input_cost = [1e10; zeros(kmax,1)];
lambda_sum = [1e10; zeros(kmax,1)];
total_cost = [1e20; zeros(kmax,1)];

% output
quiet = 0;

% Set defaults for cvx
cvx_solver gurobi
cvx_precision default


% input initilizations

U_a_p = zeros(3*time_horizon,1);
U_b_p = zeros(3*time_horizon,1);
U_c_p = zeros(3*time_horizon,1);
U_d_p = zeros(3*time_horizon,1);
U_e_p = zeros(3*time_horizon,1);

mean_X_a_no_input = Ad_concat * x_0_a;
mean_X_b_no_input = Ad_concat * x_0_b;
mean_X_c_no_input = Ad_concat * x_0_c;
mean_X_d_no_input = Ad_concat * x_0_d;
mean_X_e_no_input = Ad_concat * x_0_e;
mean_X_a = mean_X_a_no_input;
mean_X_b = mean_X_b_no_input;
mean_X_c = mean_X_c_no_input;
mean_X_d = mean_X_d_no_input;
mean_X_e = mean_X_e_no_input;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = tic;
while k <= kmax 

    % update collision avoid probabilities and gradient
    [g_ab, del_g_ab] = update_g(mean_X_a, mean_X_b, Bd_concat, time_horizon-1, r, 6);
    [g_ac, del_g_ac] = update_g(mean_X_a, mean_X_c, Bd_concat, time_horizon-1, r, 6);
    [g_ad, del_g_ad] = update_g(mean_X_a, mean_X_d, Bd_concat, time_horizon-1, r, 6);
    [g_ae, del_g_ae] = update_g(mean_X_a, mean_X_e, Bd_concat, time_horizon-1, r, 6);
    [g_bc, del_g_bc] = update_g(mean_X_b, mean_X_c, Bd_concat, time_horizon-1, r, 6);
    [g_bd, del_g_bd] = update_g(mean_X_b, mean_X_d, Bd_concat, time_horizon-1, r, 6);
    [g_be, del_g_be] = update_g(mean_X_b, mean_X_e, Bd_concat, time_horizon-1, r, 6);
    [g_cd, del_g_cd] = update_g(mean_X_c, mean_X_d, Bd_concat, time_horizon-1, r, 6);
    [g_ce, del_g_ce] = update_g(mean_X_c, mean_X_e, Bd_concat, time_horizon-1, r, 6);
    [g_de, del_g_de] = update_g(mean_X_d, mean_X_e, Bd_concat, time_horizon-1, r, 6);

    cvx_begin quiet
        variable U_a(3 * time_horizon,1);
        variable U_b(3 * time_horizon,1);
        variable U_c(3 * time_horizon,1);
        variable U_d(3 * time_horizon,1);
        variable U_e(3 * time_horizon,1);

        variable mean_X_a(6 * time_horizon, 1);
        variable mean_X_b(6 * time_horizon, 1);
        variable mean_X_c(6 * time_horizon, 1);
        variable mean_X_d(6 * time_horizon, 1);
        variable mean_X_e(6 * time_horizon, 1);

        variable lambda_ab(time_horizon-1, 1);
        variable lambda_ac(time_horizon-1, 1);
        variable lambda_ad(time_horizon-1, 1);
        variable lambda_ae(time_horizon-1, 1);
        variable lambda_bc(time_horizon-1, 1);
        variable lambda_bd(time_horizon-1, 1);
        variable lambda_be(time_horizon-1, 1);
        variable lambda_cd(time_horizon-1, 1);
        variable lambda_ce(time_horizon-1, 1);
        variable lambda_de(time_horizon-1, 1);

        variable gamma_ab(time_horizon-1, 1);
        variable gamma_ac(time_horizon-1, 1);
        variable gamma_ad(time_horizon-1, 1);
        variable gamma_ae(time_horizon-1, 1);
        variable gamma_bc(time_horizon-1, 1);
        variable gamma_bd(time_horizon-1, 1);
        variable gamma_be(time_horizon-1, 1);
        variable gamma_cd(time_horizon-1, 1);
        variable gamma_ce(time_horizon-1, 1);
        variable gamma_de(time_horizon-1, 1);
        
        variable beta_approx_ab(time_horizon-1, 1);
        variable beta_approx_ac(time_horizon-1, 1);
        variable beta_approx_ad(time_horizon-1, 1);
        variable beta_approx_ae(time_horizon-1, 1);
        variable beta_approx_bc(time_horizon-1, 1);
        variable beta_approx_bd(time_horizon-1, 1);
        variable beta_approx_be(time_horizon-1, 1);
        variable beta_approx_cd(time_horizon-1, 1);
        variable beta_approx_ce(time_horizon-1, 1);
        variable beta_approx_de(time_horizon-1, 1);

        variable delta_a(n_lin_state_a, 1);
        variable delta_b(n_lin_state_b, 1);
        variable delta_c(n_lin_state_b, 1);
        variable delta_d(n_lin_state_b, 1);
        variable delta_e(n_lin_state_c, 1);

        variable t_approx_a(n_lin_state_a, 1);
        variable t_approx_b(n_lin_state_b, 1);
        variable t_approx_c(n_lin_state_c, 1);
        variable t_approx_d(n_lin_state_d, 1);
        variable t_approx_e(n_lin_state_e, 1);
        
        variable sum_lambda(1,1);

        minimize (tau(k)*sum_lambda + U_a'*U_a + U_b'*U_b + U_c'*U_c + U_d'*U_d+ U_e'*U_e)
        subject to
            %----------------------------
            % cost variables
            %----------------------------
            sum_lambda == sum(lambda_ab) + sum(lambda_ac) + sum(lambda_ad) + sum(lambda_ae) + ...
                                           sum(lambda_bc) + sum(lambda_bd) + sum(lambda_be) + ...
                                                            sum(lambda_cd) + sum(lambda_ce) + ...
                                                                             sum(lambda_de);
            
            %----------------------------
            % linear equations defining the state
            %----------------------------
            mean_X_a == mean_X_a_no_input + Bd_concat * U_a;
            mean_X_b == mean_X_b_no_input + Bd_concat * U_b;
            mean_X_c == mean_X_c_no_input + Bd_concat * U_c;
            mean_X_d == mean_X_d_no_input + Bd_concat * U_d;
            mean_X_e == mean_X_e_no_input + Bd_concat * U_e; 

            %----------------------------
            % u \in \mathcal{U} 
            %----------------------------
            input_space_A * U_a <= input_space_b;
            input_space_A * U_b <= input_space_b; 
            input_space_A * U_c <= input_space_b;
            input_space_A * U_d <= input_space_b;
            input_space_A * U_e <= input_space_b;

            %----------------------------
            % colission avoidance constraint
            %----------------------------

            % difference of convex function representation of 
            % ||x_a - x_b||^2 >= (r + \hat{chi}*min_eig(Sigma_w(k))^1/2)^2 - slack
            % slack variables added for feasibility.
            lambda_ab >= 0;
            lambda_ac >= 0;
            lambda_ad >= 0;
            lambda_ae >= 0;
            lambda_bc >= 0;
            lambda_bd >= 0;
            lambda_be >= 0;
            lambda_cd >= 0;
            lambda_ce >= 0;
            lambda_de >= 0;
            
            
            for gamma_indx = 1:(time_horizon-1)
                beta_approx_ab(gamma_indx) >= beta_invcdf_m.* gamma_ab(gamma_indx) + beta_invcdf_c;
                beta_approx_ac(gamma_indx) >= beta_invcdf_m.* gamma_ac(gamma_indx) + beta_invcdf_c;
                beta_approx_ad(gamma_indx) >= beta_invcdf_m.* gamma_bc(gamma_indx) + beta_invcdf_c;
                beta_approx_ae(gamma_indx) >= beta_invcdf_m.* gamma_ab(gamma_indx) + beta_invcdf_c;
                beta_approx_bc(gamma_indx) >= beta_invcdf_m.* gamma_ac(gamma_indx) + beta_invcdf_c;
                beta_approx_bd(gamma_indx) >= beta_invcdf_m.* gamma_bc(gamma_indx) + beta_invcdf_c;
                beta_approx_be(gamma_indx) >= beta_invcdf_m.* gamma_ab(gamma_indx) + beta_invcdf_c;
                beta_approx_cd(gamma_indx) >= beta_invcdf_m.* gamma_ac(gamma_indx) + beta_invcdf_c;
                beta_approx_ce(gamma_indx) >= beta_invcdf_m.* gamma_bc(gamma_indx) + beta_invcdf_c;
                beta_approx_de(gamma_indx) >= beta_invcdf_m.* gamma_ab(gamma_indx) + beta_invcdf_c;
            end

            gamma_ab >= lb_approx;
            gamma_ac >= lb_approx;
            gamma_ad >= lb_approx;
            gamma_ae >= lb_approx;
            gamma_bc >= lb_approx;
            gamma_bd >= lb_approx;
            gamma_be >= lb_approx;
            gamma_cd >= lb_approx;
            gamma_ce >= lb_approx;
            gamma_de >= lb_approx;

            gamma_ab <= ub_approx;
            gamma_ac <= ub_approx;
            gamma_ad <= ub_approx;
            gamma_ae <= ub_approx;
            gamma_bc <= ub_approx;
            gamma_bd <= ub_approx;
            gamma_be <= ub_approx;
            gamma_cd <= ub_approx;
            gamma_ce <= ub_approx;
            gamma_de <= ub_approx;
            
            g_ab + del_g_ab * [U_a - U_a_p;U_b - U_b_p] - beta_approx_ab .* max_sigma + lambda_ab >= ...
                 0;
            g_ac + del_g_ac * [U_a - U_a_p;U_c - U_c_p] - beta_approx_ac .* max_sigma + lambda_ac >= ...
                 0;
            g_ad + del_g_ad * [U_a - U_a_p;U_d - U_d_p] - beta_approx_ad .* max_sigma + lambda_ad >= ...
                 0;
            g_ae + del_g_ae * [U_a - U_a_p;U_e - U_e_p] - beta_approx_ae .* max_sigma + lambda_ae >= ...
                 0;
            g_bc + del_g_bc * [U_b - U_b_p;U_c - U_c_p] - beta_approx_bc .* max_sigma + lambda_bc >= ...
                 0;
            g_bd + del_g_bd * [U_b - U_b_p;U_d - U_d_p] - beta_approx_bd .* max_sigma + lambda_bd >= ...
                 0;
            g_be + del_g_be * [U_b - U_b_p;U_e - U_e_p] - beta_approx_be .* max_sigma + lambda_be >= ...
                 0;
            g_cd + del_g_cd * [U_c - U_c_p;U_d - U_d_p] - beta_approx_cd .* max_sigma + lambda_cd >= ...
                 0;
            g_ce + del_g_ce * [U_c - U_c_p;U_e - U_e_p] - beta_approx_ce .* max_sigma + lambda_ce >= ...
                 0;
            g_de + del_g_de * [U_d - U_d_p;U_e - U_e_p] - beta_approx_de .* max_sigma + lambda_de >= ...
                 0;
            %----------------------------
            % terminal state constraint
            %----------------------------

            % approximation of inverse normal in convex region
            for delta_indx = 1:n_lin_state_a
                t_approx_a(delta_indx) >= t_invcdf_m.* delta_a(delta_indx) + t_invcdf_c;
            end
            for delta_indx = 1:n_lin_state_b
                t_approx_b(delta_indx) >= t_invcdf_m.* delta_b(delta_indx) + t_invcdf_c;
            end
            for delta_indx = 1:n_lin_state_c
                t_approx_c(delta_indx) >= t_invcdf_m.* delta_c(delta_indx) + t_invcdf_c;
            end
            for delta_indx = 1:n_lin_state_d
                t_approx_d(delta_indx) >= t_invcdf_m.* delta_d(delta_indx) + t_invcdf_c;
            end
            for delta_indx = 1:n_lin_state_e
                t_approx_e(delta_indx) >= t_invcdf_m.* delta_e(delta_indx) + t_invcdf_c;
            end

            % \mu_v in target shrunk by \beta
            target_set_a.A * mean_X_a(end-5:end) + scaled_sigma_a_vec .* t_approx_a - target_set_a.b <= 0;
            target_set_b.A * mean_X_b(end-5:end) + scaled_sigma_b_vec .* t_approx_b - target_set_b.b <= 0;
            target_set_c.A * mean_X_c(end-5:end) + scaled_sigma_c_vec .* t_approx_c - target_set_c.b <= 0;
            target_set_d.A * mean_X_d(end-5:end) + scaled_sigma_d_vec .* t_approx_c - target_set_d.b <= 0;
            target_set_e.A * mean_X_e(end-5:end) + scaled_sigma_e_vec .* t_approx_c - target_set_e.b <= 0;

            % \delta_i,v not infinity
            delta_a >= lb_approx;
            delta_b >= lb_approx;
            delta_c >= lb_approx;
            delta_d >= lb_approx;
            delta_e >= lb_approx;

            % \delta_i,v in convex region
            delta_a <= ub_approx;
            delta_b <= ub_approx;
            delta_c <= ub_approx;
            delta_d <= ub_approx;
            delta_e <= ub_approx;

            %----------------------------
            % overall safety
            %----------------------------
            sum(delta_a) + sum(delta_b) + sum(delta_c) + sum(delta_d) + sum(delta_e) <= 1 - safety_target;
            sum(gamma_ab) + sum(gamma_ac) + sum(gamma_ad) + sum(gamma_ae) + ...
                            sum(gamma_bc) + sum(gamma_bd) + sum(gamma_be) + ...
                                            sum(gamma_cd) + sum(gamma_ce) + ...
                                                            sum(gamma_de) <= 1 - safety_collision;
    cvx_end

    % update Costs
    input_cost(k+1) = U_a'*U_a + U_b'*U_b + U_c'*U_c + U_d'*U_d+ U_e'*U_e;
    lambda_sum(k+1) = sum_lambda;
    total_cost(k+1) = cvx_optval;

    % calculate convergence criteria
    conv_check = abs(input_cost(k+1) - input_cost(k) + tau(k)*(lambda_sum(k+1) - lambda_sum(k)));

    % print statistics
    if ~ quiet
        fprintf('iteration: %d ', k);
        fprintf('\t %f', cvx_optval);
        fprintf('\t %e', conv_check); 
        fprintf('\t %e', lambda_sum(k+1));
        fprintf('\t %s', cvx_status);
        fprintf('\t %s \n', toc(start_time));
    end

    % check for solved status
    if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
        % check for convergence
        if (conv_check <= epsilon_dc) && (lambda_sum(k+1) <= epsilon_lambda)                 
           break
        end

        % if not converged update previous answer to current answer
        U_a_p = U_a;
        U_b_p = U_b;
        U_c_p = U_c;
        U_d_p = U_d;
        U_e_p = U_e;

    % if results are NaN break before error
    elseif strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
        break
    end

    % update itteration number
     k = k + 1;
end
total_time = toc(start_time);

% make k not less than or equal to kmax
k = min(k, kmax);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% print some useful information
%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s \n', cvx_status);
fprintf('Computation time (sec): %f \n', total_time);
fprintf('Total Cost: %f \n', total_cost(k+1));
fprintf('Slack Cost: %f \n', lambda_sum(k+1));
fprintf('Input Cost: %f \n', input_cost(k+1));



%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make some plots
%%%%%%%%%%%%%%%%%%%%%%%%%%

F = [1, 2, 3, 4;
     1, 2, 7, 8;
     1, 4, 6, 8;
     2, 3, 5, 7;
     3, 4, 6, 5;
     5, 6, 8, 7];

red = [224 0 0] ./ 255;
blue = [30 144 255] ./ 255;
green = [0 170 85] ./ 255;
purple = [118 0 168] ./ 255;
grey = [46 52 59] ./ 255;

fig = figure();
fig.Units    = 'inches';
fig.Position = [0.75,-1,13.5,11.5];
hold on
pa = plot3([x_0_a(1); mean_X_a(1:6:end)], [x_0_a(2); mean_X_a(2:6:end)], [x_0_a(3); mean_X_a(3:6:end)], 'Color', red, 'Marker', 'h', 'LineWidth', 1);
pb = plot3([x_0_b(1); mean_X_b(1:6:end)], [x_0_b(2); mean_X_b(2:6:end)], [x_0_b(3); mean_X_b(3:6:end)], 'Color', blue, 'Marker', 'p', 'LineWidth', 1);
pc = plot3([x_0_c(1); mean_X_c(1:6:end)], [x_0_c(2); mean_X_c(2:6:end)], [x_0_c(3); mean_X_c(3:6:end)], 'Color', green, 'Marker', '^', 'LineWidth', 1);
pd = plot3([x_0_d(1); mean_X_d(1:6:end)], [x_0_d(2); mean_X_d(2:6:end)], [x_0_d(3); mean_X_d(3:6:end)], 'Color', purple, 'Marker', 'o', 'LineWidth', 1);
pe = plot3([x_0_e(1); mean_X_e(1:6:end)], [x_0_e(2); mean_X_e(2:6:end)], [x_0_e(3); mean_X_e(3:6:end)], 'Color', grey, 'Marker', 's', 'LineWidth', 1);
p1 = patch('Faces',F,'Vertices', Polyhedron(target_set_a.A([1;2;3;7;8;9],1:3), target_set_a.b([1;2;3;7;8;9])).V,...
    'FaceColor', red, ...
    'FaceAlpha',0.1); 
patch('Faces',F,'Vertices', Polyhedron(target_set_b.A([1;2;3;7;8;9],1:3), target_set_b.b([1;2;3;7;8;9])).V,...
    'FaceColor', blue, ...
    'FaceAlpha',0.1); 
patch('Faces',F,'Vertices', Polyhedron(target_set_c.A([1;2;3;7;8;9],1:3), target_set_c.b([1;2;3;7;8;9])).V,...
    'FaceColor', green, ...
    'FaceAlpha',0.1); 
patch('Faces',F,'Vertices', Polyhedron(target_set_d.A([1;2;3;7;8;9],1:3), target_set_d.b([1;2;3;7;8;9])).V,...
    'FaceColor', purple, ...
    'FaceAlpha',0.1); 
patch('Faces',F,'Vertices', Polyhedron(target_set_e.A([1;2;3;7;8;9],1:3), target_set_e.b([1;2;3;7;8;9])).V,...
    'FaceColor', grey, ...
    'FaceAlpha',0.1); 
p2 = plot3(x_0_a(1), x_0_a(2), x_0_a(3), 'Color', 'k', 'Marker', 'h', 'MarkerFaceColor', 'k', 'LineStyle','none');
plot3(x_0_b(1), x_0_b(2), x_0_b(3), 'Color', 'k', 'Marker', 'p', 'MarkerFaceColor', 'k');
plot3(x_0_c(1), x_0_c(2), x_0_c(3), 'Color', 'k', 'Marker', '^', 'MarkerFaceColor', 'k');
plot3(x_0_d(1), x_0_d(2), x_0_d(3), 'Color', 'k', 'Marker', 'o', 'MarkerFaceColor', 'k');
plot3(x_0_e(1), x_0_e(2), x_0_e(3), 'Color', 'k', 'Marker', 's', 'MarkerFaceColor', 'k');

xlabel('x (in meters)')
ylabel('y (in meters)')
zlabel('z (in meters)')
drawnow()
axis([-20 120 -35 35 -35 35])
hold off
set(gca, 'OuterPosition', [0.025,(1-0.95/2)/2,0.95,0.95/2]);

view(3)

l = legend([pa,pb,pc,pd,pe,p1,p2], {'A', 'B', 'C', 'D', 'E', 'Target Set', 'Initial Location' },...
    'Orientation','horizontal', ...
    'Location', 'northoutside', ...
    'NumColumns', 4, ...
    'interpreter', 'latex');
