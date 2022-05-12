function [P_target, P_collision] = verify_relative(mean_traj, Wd_concat, psi, nu, time_horizon, target_sets, target_times, r, samples)
        
    
    rng('default');
    disturbed_mean = mvnrnd(zeros(size(psi,1),1),psi,samples)'./sqrt(chi2rnd(nu,1,samples)./nu);
   
    
    for i = 1:samples
        disturbed_mean(:,i) = mean_traj +  Wd_concat * disturbed_mean(:,i);
    end
   
    targets = max(size(target_times));
    in_target = zeros(targets, samples);
	for i = 1:targets
        index = 6*(target_times(i)-1)+[1:6];
        for j = 1:samples
            in_target(i, j) = target_sets(i).contains( disturbed_mean(index, j) );
        end
    end
    P_target = sum(sum(in_target, 1) == targets) / samples;
    
    
    collision_1v = zeros(time_horizon, samples);
    for t = 1:time_horizon
        time_index = 6*(t-1) + [1:2];
        for k = 1:samples
                collision_1v(t, k) = ( norm( disturbed_mean(time_index, k)) >= r); 
        end
    end
    P_collision = sum(sum(collision_1v, 1)  == time_horizon) / samples; 
end