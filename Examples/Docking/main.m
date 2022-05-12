%% clean environment
clear
clc
close all
cvx_clear;

addpath '../'

%% setup system

sampling_period              = 60*5;                                            % sec
orbital_radius               = (35622 + 6378.1) * 1000;                         % m
gravitational_constant       = 6.673e-11;                                       % m^3 kg^-1 sec^-2
celestial_mass               = 5.9472e24;                                       % kg
gravitational_body           = gravitational_constant * celestial_mass;         % m^3 sec^-2
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
%%
% Impulse control
Bd = Ad*[zeros(3); eye(3)];

%% problem set up
time_horizon = 8;

% initial states
% format: x, y, z,  x., y., z.
x_0_a = [90;  -5;  0.1; 0; 0; 0] ; % satellite A
x_0_b = [95;   5; -0.1; 0; 0; 0] ; % satellite B
x_0_c = [100; -5; -0.1; 0; 0; 0] ; % satellite C
x_0_d = [105;  5; -0.1; 0; 0; 0] ; % satellite E
x_0_e = [110; -5;  0.1; 0; 0; 0] ; % satellite D
x_0_f = [ 95;  0;   10; 0; 0; 0] ; % satellite F
x_0_g = [ 95;  0;  -10; 0; 0; 0] ; % satellite G
x_0 = [x_0_a, x_0_b, x_0_c, x_0_d, x_0_e, x_0_f, x_0_g];
vehicles = 7;
combinations = vehicles * (vehicles-1) / 2;

% target sets
% format: x, y, z, x., y., z.

target_set_a = Polyhedron('lb', [ -2.5; -2.5;  7.5; -0.01; -0.01; -0.01], ...
                          'ub', [  2.5;  2.5; 12.5;  0.01;  0.01;  0.01]);  
target_set_b = Polyhedron('lb', [ -2.5; -2.5;-12.5; -0.01; -0.01; -0.01], ...
                          'ub', [  2.5;  2.5; -7.5;  0.01;  0.01;  0.01]);  
target_set_c = Polyhedron('lb', [-7.5;     5; -2.5; -0.01; -0.01; -0.01], ... 
                          'ub', [-2.5;    10;  2.5;  0.01; 0.01;  0.01]);   
target_set_d = Polyhedron('lb', [-7.5;   -10; -2.5; -0.01; -0.01; -0.01], ... 
                          'ub', [-2.5;    -5;  2.5;  0.01; 0.01;  0.01]);   
target_set_e = Polyhedron('lb', [ 2.5;   2.5; -2.5; -0.01; -0.01; -0.01], ... 
                          'ub', [ 7.5;   7.5;  2.5;  0.01; 0.01;  0.01]);   
target_set_f = Polyhedron('lb', [ 2.5;  -7.5; -2.5; -0.01; -0.01; -0.01], ...
                          'ub', [ 7.5;  -2.5;  2.5;  0.01; 0.01;  0.01]);   
target_set_g = Polyhedron('lb', [-12.5; -2.5; -2.5; -0.01; -0.01; -0.01], ...
                          'ub', [-7.5;   2.5;  2.5;  0.01; 0.01;  0.01]);   
                      
target_set_all_A = target_set_a.A;
target_set_all_b = [target_set_a.b, target_set_b.b, target_set_c.b, target_set_d.b, target_set_e.b, target_set_f.b, target_set_g.b];
                      
% get number of half space constraints                      
n_lin_state = size(target_set_all_A,1);
                      
% Input space
u_max = 3;
input_space = Polyhedron('lb', [-u_max; -u_max; -u_max], ... 
                         'ub', [ u_max;  u_max;  u_max]);                         

input_space_A = blkdiag(input_space.A);
for i=1:(time_horizon-1)
    input_space_A = blkdiag(input_space_A, input_space.A);
end

input_space_b = repmat(input_space.b, time_horizon,1);

% collision avoid region radius
r = 8;

% safety threshold
safety_target           = 0.8;  % in target set
safety_collision_2_v    = 0.8;  % intersatellite
safety_collision_1_v    = 0.8;  % between satellite and chief



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
nu = 20;

Cw_psi_Cw = Wd_concat * psi * Wd_concat';

% multiplier for target set (same for all vehicles)
scaled_sigma_vec = diag(target_set_a.A * Cw_psi_Cw(end-5:end,end-5:end) * target_set_a.A').^(1/2);

% multiplier for collision avoidance
max_sigma = zeros(time_horizon-1,1);
for i = 1:(time_horizon-1)
    index = 6*(i-1) + (1:3); 
    max_sigma(i) = sqrt(nu * max(eig(Cw_psi_Cw(index, index))));
end

% disturbance  linearization
%set up
syms x t_pdf(x) beta_prime_pdf(x, a, b);

by = 5e-6; % step size
approx_to = 0.99999;
approx_from = .5;
tolerance = 1e-2;
linearize_from = safety_collision_2_v;

lb_approx = 1 - approx_to;
ub_approx = 1 - linearize_from;

% t dist
t_init = 0; % known quantile at p=0.5
t_pdf(x) = gamma((nu+1)/2) / (sqrt(nu*pi) * gamma(nu/2)) * (1 + x^2/nu)^(-(nu+1)/2);
[t_invcdf_m, t_invcdf_c] = quantile_approx(1, by, approx_to, approx_from, t_init, t_pdf, tolerance, linearize_from);


% beta prime setup
med_approx_A = log(2)-1/3;
med_approx_B = 1;
beta_prime_med_approx = @(a,b) sqrt((2^(-1/a)*( med_approx_A + med_approx_B*a))  / (2^(-1/b)*( med_approx_A + med_approx_B*b)));
beta_prime_pdf(x, a, b) =  2*x^(2*a-1) * (1+x^2)^(-a-b) / beta(a,b);

% beta prime dist collision avoid with chief
q = size(Ad, 1)/2;
a_1v =  q / 2;
b_1v = nu / 2;
[beta_invcdf_m_1v, beta_invcdf_c_1v] = quantile_approx(1, by, approx_to, approx_from, beta_prime_med_approx(a_1v,b_1v), beta_prime_pdf(x, a_1v, b_1v), tolerance, linearize_from);

% beta prime dist intervehicle collision avoid
a_2v = 2 * q * ( q/2 + nu^2/4 - nu + q*nu/2 - 2*q + 1) / ((nu - 2)*(q + nu - 2));
b_2v = 2 * (-q + nu^2/4 - nu/2 + q*nu/2) / (q + nu - 2);
[beta_invcdf_m_2v, beta_invcdf_c_2v] = quantile_approx(1, by, approx_to, approx_from, beta_prime_med_approx(a_2v,b_2v), beta_prime_pdf(x, a_2v, b_2v), tolerance, linearize_from);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
kmax = 100;

% convergence perameters
epsilon_dc = 1e-8; % convergence in cost
epsilon_lambda = 1e-8; % convergence of sum of slack variables to zero

% cost of slack variable
tau_max = 1000;
gamma = 1.2;
tau = min(tau_max * ones(kmax,1),  gamma.^(0:(kmax-1))');

% storage initial cost for convergence check
input_cost = [1e10; zeros(kmax,1)];
lambda_sum = [1e10; zeros(kmax,1)];
total_cost = [1e20; zeros(kmax,1)];

% output
quiet = 0;

% Set defaults for cvx
cvx_solver gurobi
cvx_precision default

% initial input guess
U_p = zeros(3*time_horizon, vehicles);

% state propigated without input
mean_X_no_input = Ad_concat * x_0;
mean_X = mean_X_no_input + Bd_concat * U_p;


%holders for norm gradient approximation
S = [eye(3), zeros(3)];

norm_approx_2v = zeros(time_horizon-1, combinations);
norm_approx_gradient_2v = zeros(time_horizon-1, 2*size(Bd_concat,2), combinations);

norm_approx_1v = zeros(time_horizon-1, vehicles);
norm_approx_gradient_1v = zeros(time_horizon-1, size(Bd_concat,2), vehicles);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
while k <= kmax 

    % update collision avoid gradient
    
    for i = 1:(vehicles-1)
        for j = (i+1):vehicles
            index = (i-1)*(vehicles-1-i/2) + j-1;
            [norm_approx_2v(:,index), norm_approx_gradient_2v(:,:,index)] = ...
                update_g_2v(mean_X(:,i), mean_X(:,j), Bd_concat, time_horizon-1, r, S);
        end
    end
    
    for i = 1:vehicles
        [norm_approx_1v(:,i), norm_approx_gradient_1v(:,:,i)] = ...
            update_g_1v(mean_X(:,i), Bd_concat, time_horizon-1, r, S);
    end
    
    cvx_begin quiet
        variable U(3 * time_horizon, vehicles);
        variable mean_X(6 * time_horizon, vehicles);

        % 2 vehicle collision avoid
        variable lambda_2v(time_horizon-1, combinations) nonnegative;
        variable gamma_2v(time_horizon-1, combinations);
        variable beta_approx_2v(time_horizon-1, combinations);

        % collision avoid with chief
        variable lambda_1v(time_horizon-1, vehicles) nonnegative;
        variable gamma_1v(time_horizon-1, vehicles);
        variable beta_approx_1v(time_horizon-1, vehicles);
        
        % targe set
        variable delta(n_lin_state, vehicles);
        variable t_approx(n_lin_state, vehicles);
        
        % cost params
        variable sum_lambda(1,1);
        variable quad_input_cost(1,1);

        minimize (tau(k)*sum_lambda + quad_input_cost)
        subject to
            %----------------------------
            % cost variables
            %----------------------------
            sum_lambda == sum(vec(lambda_2v)) + sum(sum(lambda_1v,1),2); 
                      
            quad_input_cost >=  vec(U)'*vec(U);
            
            %----------------------------
            % linear equations defining the state
            %----------------------------
            mean_X == mean_X_no_input + Bd_concat * U;

            %----------------------------
            % u \in \mathcal{U} 
            %----------------------------
            for i = 1:vehicles
                input_space_A * U(:,i) <= input_space_b;
            end
            
            %----------------------------
            % colission avoidance constraint (intervehicle)
            %----------------------------
            
            % quantile in region of interest
            vec(gamma_2v) >= lb_approx;
            vec(gamma_2v) <= ub_approx;

            for i = 1:(vehicles-1)
                for j = (i+1):vehicles
                    index = (i-1)*(vehicles-1-i/2) + j-1;
                                
                    % quantile approx
                    for gamma_indx = 1:(time_horizon-1)
                        beta_approx_2v(gamma_indx, index) >= beta_invcdf_m_2v.* gamma_2v(gamma_indx, index) + beta_invcdf_c_2v;
                    end
                                                                                                                
                    % difference of convex collision avoidance constraints
                    norm_approx_2v(:,index) + ...
                        norm_approx_gradient_2v(:,:,index) * [U(:,i) - U_p(:,i);U(:,j) - U_p(:,j)] - ...
                        beta_approx_2v(:,index) .* (sqrt(2)*max_sigma) + lambda_2v(:,index) >= 0;

                end
            end
            %----------------------------
            % colission avoidance constraint (from chief)
            %----------------------------
                                
            % quantile in region of interest
            vec(gamma_1v) >= lb_approx;
            vec(gamma_1v) <= ub_approx;
            
            for i = 1:vehicles
            
                    % quantile approx
                    for gamma_indx = 1:(time_horizon-1)
                        beta_approx_1v(gamma_indx, i) >= beta_invcdf_m_1v.* gamma_1v(gamma_indx, i) + beta_invcdf_c_1v;
                    end
                    
                    % difference of convex collision avoidance constraints
                    norm_approx_1v(:,i) + ...
                        norm_approx_gradient_1v(:,:,i) * (U(:,i) - U_p(:,i)) - ...
                        beta_approx_1v(:,i) .* max_sigma + lambda_1v(:,i) >= 0;
            end
            
            %----------------------------
            % terminal state constraint
            %----------------------------
            % quantile in reqion of interest
            vec(delta) >= lb_approx;
            vec(delta) <= ub_approx;
            
            for i = 1:vehicles
                % quantile approx
                for delta_indx = 1:n_lin_state
                    t_approx(delta_indx, i) >= t_invcdf_m.* delta(delta_indx, i) + t_invcdf_c;
                end
                
                % mean in shrunk target set
                target_set_all_A * mean_X(end-5:end, i) + scaled_sigma_vec .* t_approx(:,i) - target_set_all_b(:,i) <= 0;
            end     
            %----------------------------
            % overall safety
            %----------------------------
            sum(vec(delta))  <= 1 - safety_target;
            sum(vec(gamma_2v)) <= 1 - safety_collision_2_v;
            sum(vec(gamma_1v)) <= 1 - safety_collision_1_v;
    cvx_end

    % update Costs
    input_cost(k+1) = quad_input_cost;
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
        fprintf('\t %f \n', toc(start_time));
    end

    % check for solved status
    if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
        % check for convergence
        if (conv_check <= epsilon_dc) && (lambda_sum(k+1) <= epsilon_lambda)                 
           break
        end

        % if not converged update previous answer to current answer
        U_p = U;

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

if strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%

target_sets(1) = target_set_a;
target_sets(2) = target_set_b;
target_sets(3) = target_set_c;
target_sets(4) = target_set_d;
target_sets(5) = target_set_e;
target_sets(6) = target_set_f;
target_sets(7) = target_set_g;
make_relative_plots

%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%

[P_target, P_2v, P_1v] = verify_relative(mean_X, Wd_concat, psi, nu, time_horizon-1, target_sets, r, 10000)
