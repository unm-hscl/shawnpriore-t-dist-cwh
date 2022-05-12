%%%%%%%%%%%%%%%%%%%%%%%%
% Blackmore 2011 for goal achievement and collision avoidance
% Primary Coder: Vignesh Sivaramakrishnan
% Modified by: Shawn Priore
%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
% big M arbitrary constant
large_constant = 5000;
N = 25;

% polytoupe defining ||x_i|| = r
Avoid_A = [  eye(2);
            -eye(2);
            1,1;
            1,-1;
            -1,1;
            -1,-1];
Avoid_b = [r * ones(4,1); r * sqrt(2) * ones(4,1)];


% randomly generate the disturbance vector from the multivariate t.
W_deputy = mvnrnd(zeros(size(psi,1),1),psi,N)'./sqrt(chi2rnd(nu,1,N)./nu);

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T.

tic;
cvx_begin 
    variable U_pc(3 * time_horizon,1);
    variable mean_X_pc(6 * time_horizon, 1);
    % incorporate disturbances 
    variable x_pc(6 * time_horizon,N);
    % target set
    variable t_pc(N) binary;
    % collision avoidance between vehicle and chief
    variable c_pc(time_horizon * size(Avoid_A,1), N) binary;
    variable sc_pc(time_horizon , N) binary;
    variable ssc_pc(N) binary;

    minimize (U_pc'*U_pc);

    subject to
        mean_X_pc == mean_X_no_input + Bd_concat * U_pc;
        x_pc(:,1:N) == Wd_concat * W_deputy + repmat(mean_X_pc,1,N);

        input_space_A * U_pc <= input_space_b;

        for i = 1:N

            concat_target_A * x_pc(:,i) - concat_target_b <= large_constant*t_pc(i);
            
            for t = 1:time_horizon
                
                Avoid_A * x_pc(6*(t-1) + (1:2),i) - Avoid_b + large_constant * c_pc(size(Avoid_A,1)*(t-1) + (1:size(Avoid_A,1)) , i) >= 0; 

                sum(c_pc(size(Avoid_A,1)*(t-1) + (1:size(Avoid_A,1)) , i)) - (size(Avoid_A,1) - 1) <= sc_pc(t, i);
                
            end

            
            sum(sc_pc(:,i)) <= time_horizon * ssc_pc(i);
        end

        1/N * sum(t_pc) <= 1-safety_target;
        1/N * sum(ssc_pc) <= 1-safety_collision;
cvx_end
total_time = toc;

fprintf('%s \n', cvx_status);
fprintf('Computation time (sec): %f \n', total_time);
fprintf('Input Cost: %f \n', cvx_optval);
