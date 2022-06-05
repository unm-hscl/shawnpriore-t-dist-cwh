function [P_target, P_2v] = verify_relative(mean_traj, Wd_concat, psi, nu, time_horizon, target_sets, r, samples)
    vehicles = size(mean_traj, 2);
    combinations = vehicles * (vehicles-1) / 2;
        
    disturbed_mean = zeros([size(mean_traj),samples]);
    
    for i = 1:vehicles
        disturbed_mean(:,i,:) = mvnrnd(zeros(size(psi,1),1),psi,samples)'./sqrt(chi2rnd(nu,1,samples)./nu);
    end
    
    for i = 1:samples
        disturbed_mean(:,:,i) = mean_traj +  Wd_concat * disturbed_mean(:,:,i);
    end
   
    in_target = zeros(vehicles, samples);
    for i = 1:vehicles
        for j = 1:samples
            in_target(i, j) = target_sets(i).contains( disturbed_mean(end-5:end, i, j) );
        end
    end
    P_target = sum(sum(in_target, 1) == vehicles) / samples;
    
    collision_2v = zeros(combinations, time_horizon, samples);
    for i = 1:(vehicles-1)
        for j = (i+1):vehicles
            index = (i-1)*(vehicles-1-i/2) + j-1;
            for t = 1:time_horizon
                time_index = 6*(t-1) + [1:3];
                for k = 1:samples
                    collision_2v(index, t, k) = ( norm( disturbed_mean(time_index, i, k) - disturbed_mean(time_index, j, k)) >= r);
                end
            end
        end
    end
    P_2v = sum(sum(sum(collision_2v, 1) == combinations, 2) == time_horizon) / samples;
end