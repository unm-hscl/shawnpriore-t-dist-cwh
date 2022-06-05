function [g, del_g] = update_g_2v(mu_1, mu_2, Cu, time_horizon, r, S)
    % calculate and extract input
    mu = mu_1 - mu_2;
    
    q = size(S,1);
    dim = size(S,2);
    
    % memory holders
    g = zeros(time_horizon, 1);
    gradient_g = zeros(time_horizon, size(Cu,2));
    
    % iterate through time index
    for i = 1:time_horizon
        % get relavent indexes
        index = dim*(i-1) + (1:dim);
        
        % calculate L_2 norm of mean
        mu_i = mu(index);
        g(i) = norm(S*mu_i)-r; 
        
        % get indexed rows of controlability matrix
        Cu_i = Cu(index, :);
        
        % calculate gradient of norm
        gradient_g(i,:) =  mu_i' * (S') * S * Cu_i ./ norm(S*mu_i);
    end
    
    % compile gradient w.r.t u_1 and u_2
    del_g = [gradient_g, -gradient_g];
end

