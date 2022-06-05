
%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
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

% initial input guess
U_p = zeros(3*time_horizon, vehicles);

% state propigated without input
mean_X_no_input = Ad_concat * x_0;
mean_X = mean_X_no_input + Bd_concat * U_p;


%holders for norm gradient approximation
S = [eye(3), zeros(3)];

norm_approx_2v = zeros(time_horizon-1, combinations);
norm_approx_gradient_2v = zeros(time_horizon-1, 2*size(Bd_concat,2), combinations);

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
    
    
    cvx_begin quiet
        variable U(3 * time_horizon, vehicles);
        variable mean_X(6 * time_horizon, vehicles);

        % 2 vehicle collision avoid
        variable lambda_2v(time_horizon-1, combinations) nonnegative;
        variable gamma_2v(time_horizon-1, combinations);
        variable beta_approx_2v(time_horizon-1, combinations);
        
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
            sum_lambda == sum(vec(lambda_2v)); 
                      
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

time_holder(mc_iter) = total_time;

if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
    input_cost_holder(mc_iter) = input_cost(k+1);
    slack_cost_holder(mc_iter) = lambda_sum(k+1);
    iter_holder(mc_iter) = k;
    success_holder(mc_iter) = 1;
end
