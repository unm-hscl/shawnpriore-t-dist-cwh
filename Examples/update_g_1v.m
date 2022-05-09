function [g, del_g] = update_g_1v(mu, Cu, time_horizon, r, S)
    
    q = size(S,1);
    dim = size(S,2);

    % memory holders
    g = zeros(time_horizon, 1);
    del_g = zeros(time_horizon, size(Cu,2));
    
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
        del_g(i,:) =  mu_i' * (S') * S * Cu_i ./ norm(mu_i);
    end
end

