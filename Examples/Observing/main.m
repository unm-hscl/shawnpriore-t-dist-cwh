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

% Continuous-time LTI CWH unforced dynamics 
A_cont = [zeros(3), eye(3);
        3*orbit_ang_vel^2, 0, 0, 0, 2*orbit_ang_vel, 0;
        0, 0, 0, -2*orbit_ang_vel, 0, 0;
        zeros(1,6)];
    

% Discrete-time system is Phi(T_s) for sampling time T_s since the system is time-invariant
Ad = expm(A_cont * sampling_period);

% Impulse control
Bd = Ad*[zeros(3); eye(3)];

%% problem set up
time_horizon = 8;

% initial states
% format: x, y, z,  x., y., z.
x_0_deputy = [10;  0;  0; 0; 0; 0] ; % deputy

% target sets
% format: x, y, theta, x., y., theta.
centroid = 10/sqrt(2);
theta_range = deg2rad(20);

target_set_2  = Polyhedron('lb', [ centroid-2; -centroid-2; -3/4*pi-theta_range; -0.01; -0.01; -0.01], ...
                           'ub', [ centroid+2; -centroid+2; -3/4*pi+theta_range;  0.01;  0.01;  0.01]);  
target_set_4  = Polyhedron('lb', [ 0-2; -10-2; -pi-theta_range; -0.01; -0.01; -0.01], ...
                           'ub', [ 0+2; -10+2; -pi+theta_range;  0.01;  0.01;  0.01]);  
target_set_6  = Polyhedron('lb', [ -centroid-2; -centroid-2; -5/4*pi-theta_range; -0.01; -0.01; -0.01], ...
                           'ub', [ -centroid+2; -centroid+2; -5/4*pi+theta_range;  0.01;  0.01;  0.01]);  
target_set_8  = Polyhedron('lb', [ -10-2; 0-2; -3/2*pi-theta_range; -0.01; -0.01; -0.01], ...
                           'ub', [ -10+2; 0+2; -3/2*pi+theta_range;  0.01;  0.01;  0.01]);  
                       
target_set_A = target_set_2.A;
temp_target_A = [zeros(size(target_set_A)), target_set_A];
concat_target_A = blkdiag(temp_target_A,temp_target_A,temp_target_A,temp_target_A);
concat_target_b = [target_set_2.b; target_set_4.b; target_set_6.b; target_set_8.b];
                      
% get number of half space constraints                      
n_lin_state = size(concat_target_A,1);
                      
% Input space
u_max = 3;
torque_max = pi/2;
input_space = Polyhedron('lb', [ -u_max; -u_max; -torque_max], ... 
                          'ub', [ u_max;  u_max;  torque_max]);                         

input_space_A = blkdiag(input_space.A);
for i=1:(time_horizon-1)
    input_space_A = blkdiag(input_space_A, input_space.A);
end

input_space_b = repmat(input_space.b, time_horizon,1);

% collision avoid region radius
r = 8;

% safety threshold
safety_target       = 0.8;  % in target sets
safety_collision    = 0.8;  % between deputy and chief



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
cov = diag([10^-4, 10^-4, 10^-6, 5*10^-8, 5*10^-8, 5*10^-10]);
psi = cov;
for i=1:(time_horizon-1)
    psi = blkdiag(psi, cov);
end
nu = 4;

Cw_psi_Cw = Wd_concat * psi * Wd_concat';

% multiplier for target set
scaled_sigma_vec = diag(concat_target_A * Cw_psi_Cw * concat_target_A').^(1/2);

% multiplier for collision avoidance
max_sigma = zeros(time_horizon,1);
for i = 1:time_horizon
    index = 6*(i-1) + (1:3); 
    max_sigma(i) = sqrt(nu * max(eig(Cw_psi_Cw(index, index))));
end

% disturbance  linearization
% set up
syms t_pdf(x) beta_prime_pdf(x);

by = 5e-6; % step size
approx_to = 0.99999;
approx_from = .5;
tolerance = 1e-2;
linearize_from = min(safety_target, safety_collision);

lb_approx = 1 - approx_to;
ub_approx = 1 - linearize_from;

% t dist
t_init = 0; % known quantile at p=0.5
t_pdf(x) = gamma((nu+1)/2) / (sqrt(nu*pi) * gamma(nu/2)) * (1 + x^2/nu)^(-(nu+1)/2);
[t_invcdf_m, t_invcdf_c] = quantile_approx(1, by, approx_to, approx_from, t_init, t_pdf, tolerance, linearize_from);

% beta prime dist intervehicle collision avoid
beta_prime_a =  2 / 2;
beta_prime_b = nu / 2;

beta_prime_med = sqrt(2^(1/beta_prime_b)-1);

beta_prime_pdf(x) = 2*x^(2*beta_prime_a-1)*(1+x^2)^(-beta_prime_a-beta_prime_b) / beta(beta_prime_a, beta_prime_b);

[beta_invcdf_m, beta_invcdf_c] = quantile_approx(1, by, approx_to, approx_from, beta_prime_med, beta_prime_pdf, tolerance, linearize_from);



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
quiet = 1;

% Set defaults for cvx
cvx_solver gurobi
cvx_precision default

% initial input guess
U_p = 0.1*randn(3*time_horizon, 1);

% state propigated without input
mean_X_no_input = Ad_concat * x_0_deputy;
mean_X = mean_X_no_input + Bd_concat * U_p;

% matrix to extract position
S = [eye(2), zeros(2,4)];


%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
while k <= kmax 

    % update collision avoid gradient    
    [norm_approx, norm_approx_gradient] = update_g_1v(mean_X, Bd_concat, time_horizon, r, S);
    
    cvx_begin quiet
        variable U(3 * time_horizon, 1);
        variable mean_X(6 * time_horizon, 1);


        % collision avoid with chief
        variable lambda(time_horizon, 1) nonnegative;
        variable gamma(time_horizon, 1);
        variable beta_approx(time_horizon, 1);
        
        % targe set
        variable delta(n_lin_state, 1);
        variable t_approx(n_lin_state, 1);
        
        % cost params
        variable sum_lambda(1,1);
        variable quad_input_cost(1,1);

        minimize (tau(k)*sum_lambda + quad_input_cost)
        subject to
            %----------------------------
            % cost variables
            %----------------------------
            sum_lambda == sum(lambda);
                      
            quad_input_cost >=  U'*U;
            
            %----------------------------
            % linear equations defining the state
            %----------------------------
            mean_X == mean_X_no_input + Bd_concat * U;

            %----------------------------
            % u \in \mathcal{U} 
            %----------------------------
            input_space_A * U <= input_space_b;
                        
            %----------------------------
            % colission avoidance from chief
            %----------------------------
                                
            % quantile in region of interest
            gamma >= lb_approx;
            gamma <= ub_approx;
            
            % quantile approx
            for gamma_indx = 1:time_horizon
                beta_approx(gamma_indx) >= beta_invcdf_m.* gamma(gamma_indx) + beta_invcdf_c;
            end

            % difference of convex collision avoidance constraints
            norm_approx + ...
                norm_approx_gradient * (U - U_p) - ...
                beta_approx .* max_sigma + lambda >= 0;
            
            %----------------------------
            % terminal state constraint
            %----------------------------
            % quantile in reqion of interest
            delta >= lb_approx;
            delta <= ub_approx;
            
            % quantile approx
            for delta_indx = 1:n_lin_state
                t_approx(delta_indx) >= t_invcdf_m.* delta(delta_indx) + t_invcdf_c;
            end

            % mean in shrunk target set
            concat_target_A * mean_X + scaled_sigma_vec .* t_approx - concat_target_b <= 0;
            
            %----------------------------
            % overall safety
            %----------------------------
            sum(delta)  <= 1 - safety_target;
            sum(gamma)  <= 1 - safety_collision;
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

fprintf('\n%s in ', cvx_status);
fprintf('%i itterations \n', k);
fprintf('Computation time (sec): %f \n', total_time);
fprintf('Total Cost: %f \n', total_cost(k+1));
fprintf('Slack Cost: %f \n', lambda_sum(k+1));
fprintf('Input Cost: %f \n', input_cost(k+1));

if strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify probabilities of our method
%%%%%%%%%%%%%%%%%%%%%%%%%%

target_sets(1) = target_set_2;
target_sets(2) = target_set_4;
target_sets(3) = target_set_6;
target_sets(4) = target_set_8;
target_times = [2,4,6,8];

[P_target, P_collision] = verify_relative(mean_X, Wd_concat, psi, nu, time_horizon, target_sets, target_times, r, 1e4)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% run and verify probabilities of particle control
%%%%%%%%%%%%%%%%%%%%%%%%%%

particle_control

if strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
    return
end

[P_target_pc, P_collision_pc] = verify_relative(mean_X_pc, Wd_concat, psi, nu, time_horizon, target_sets, target_times, r, 1e4)


%%%%%%%%%%%%%%%%%%%%%%%%%%
% make nice plot
%%%%%%%%%%%%%%%%%%%%%%%%%%
make_plots
