%% clean environment
clear
clc
close all
cvx_clear;

addpath '../'

% output
quiet = 1;

% Set defaults for cvx
cvx_solver gurobi
cvx_precision high

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
x_0_not_random = [x_0_a, x_0_b, x_0_c];
vehicles = 3;
combinations = vehicles * (vehicles-1) / 2;

% target sets
% format: x, y, z, x., y., z.

target_set_a = Polyhedron('lb', [ -3; -3;  6; -0.1; -0.1; -0.1], ...
                          'ub', [  3;  3; 12;  0.1;  0.1;  0.1]);  
target_set_b = Polyhedron('lb', [ -3; -3;-12; -0.1; -0.1; -0.1], ...
                          'ub', [  3;  3; -6;  0.1;  0.1;  0.1]);  
target_set_c = Polyhedron('lb', [ -9;  6; -3; -0.1; -0.1; -0.1], ... 
                          'ub', [ -3; 12;  3;  0.1;  0.1;  0.1]);   
                      
target_set_all_A = target_set_a.A;
target_set_all_b = [target_set_a.b, target_set_b.b, target_set_c.b];
                      
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
nu = 4;

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

% beta prime dist intervehicle collision avoid
q = 3;
a_2v = 2 * q * ( q/2 + nu^2/4 - nu + q*nu/2 - 2*q + 1) / ((nu - 2)*(q + nu - 2));
b_2v = 2 * (-q + nu^2/4 - nu/2 + q*nu/2) / (q + nu - 2);
[beta_invcdf_m_2v, beta_invcdf_c_2v] = quantile_approx(1, by, approx_to, approx_from, beta_prime_med_approx(a_2v,b_2v), beta_prime_pdf(x, a_2v, b_2v), tolerance, linearize_from);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% set variable for probability verify
%%%%%%%%%%%%%%%%%%%%%%%%%%

target_sets(1) = target_set_a;
target_sets(2) = target_set_b;
target_sets(3) = target_set_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% holders
%%%%%%%%%%%%%%%%%%%%%%%%%%
mc_samples = 1000;
P_target_holder = zeros(mc_samples,1);
P_2v_holder = zeros(mc_samples,1);
time_holder = zeros(mc_samples,1);
input_cost_holder = zeros(mc_samples,1);
iter_holder = zeros(mc_samples,1);
success_holder = zeros(mc_samples,1);
x_mean_holder = zeros(6*(time_horizon+1), 3, mc_samples);

for mc_iter = 1:mc_samples
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % randomize init conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    clc
    mc_iter
    x_0 = x_0_not_random + mvnrnd(zeros(6,1),blkdiag(eye(3), zeros(3)),3)'./sqrt(chi2rnd(10,1,3)./10);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve the problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    solve
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % verify probabilities
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
        x_mean_holder(:,:,mc_iter) = [x_0; mean_X];
        [P_target_holder(mc_iter), P_2v_holder(mc_iter)] = verify_relative(mean_X, Wd_concat, psi, nu, time_horizon, target_sets, r, 10000);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_paths


descript = @(x) [mean(x), std(x), min(x), max(x)];
descript(time_holder)
descript(input_cost_holder)
descript(slack_cost_holder)
descript(iter_holder)
descript(P_2v_holder)
descript(P_target_holder)
descript(success_holder)
