colors = [224,   0,   0; % red
           30, 144  255; % dark blue
            0, 170,  85; % green
          118,   0, 168; % purple
           46,  52,  59; % grey
          236, 176,  31;  % yellow
           76, 189, 237; % light blue
          161,  19,  46  % dark red
           ] ./ 255;

plot_symbols = ['o', 's', 'd'];

       
F_xy = [1, 4, 6, 8];
F_xz = [3, 4, 6, 5];


fig = figure();
fig.WindowState = 'maximized';

subplot(7,2,[1 2]);
hold on 
for i = 1:vehicles
    plot(nan, nan, 'Color', colors(i, :), 'Marker', plot_symbols(i));
end
patch('Faces',F_xy,'Vertices', Polyhedron(target_sets(1).A([1;2;3;7;8;9],1:3), target_sets(1).b([1;2;3;7;8;9])).V,...
    'FaceColor',  'white', ...
    'EdgeColor', 'black', ...
    'FaceAlpha', 0); 
    errorbar(10,10,1,1,1,1,'--',  'Color', 'k')
plot(x_0(1,1), nan, 'Color', 'k', 'Marker', plot_symbols(1), 'MarkerFaceColor', 'k', 'LineStyle','none');
plots=get(gca, 'Children');

legend([plots(6), plots(5), plots(4), plots(2), plots(3), plots(1)], ...
     {'A', 'B', 'C', 'Trajectory Bounds', 'Target Set', 'Non Random Initial Location' },...
    'Orientation','horizontal', ...
    'Location', 'south', ...
    'NumColumns', 6, ...
    'interpreter', 'latex');

axis([0 1 0 1]);
axis off
hold off

subplot(7,2,[3:2:14]);
hold on
for veh = 1:vehicles
    x_mean = mean(reshape(x_mean_holder(1:6:end,veh,:), time_horizon+1, mc_samples),2);
    y_mean =  mean(reshape(x_mean_holder(2:6:end,veh,:), time_horizon+1, mc_samples),2);
    xneg = min(reshape(x_mean_holder(1:6:end,veh,:), time_horizon+1, mc_samples),[],2)-x_mean;
    xpos = max(reshape(x_mean_holder(1:6:end,veh,:), time_horizon+1, mc_samples),[],2)-x_mean;
    yneg = min(reshape(x_mean_holder(2:6:end,veh,:), time_horizon+1, mc_samples),[],2)-y_mean;
    ypos = max(reshape(x_mean_holder(2:6:end,veh,:), time_horizon+1, mc_samples),[],2)-y_mean;
    errorbar(x_mean,y_mean,yneg,ypos,xneg,xpos,'--',  'Color', colors(veh,:))
    patch('Faces',F_xy,'Vertices', Polyhedron(target_sets(veh).A([1;2;3;7;8;9],1:3), target_sets(veh).b([1;2;3;7;8;9])).V,...
        'FaceColor', colors(veh, :), ...
        'FaceAlpha', 0.1); 
    plot(x_0_not_random(1,veh), x_0_not_random(2,veh), 'Color', 'k', 'Marker', plot_symbols(veh), 'MarkerFaceColor', 'k', 'LineStyle','none');
end
xlabel('x (in meters)')
ylabel('y (in meters)')
hold off

subplot(7,2,[4:2:14]);
hold on
for veh = 1:vehicles
    x_mean = mean(reshape(x_mean_holder(1:6:end,veh,:), time_horizon+1, mc_samples),2);
    y_mean =  mean(reshape(x_mean_holder(3:6:end,veh,:), time_horizon+1, mc_samples),2);
    xneg = min(reshape(x_mean_holder(1:6:end,veh,:), time_horizon+1, mc_samples),[],2)-x_mean;
    xpos = max(reshape(x_mean_holder(1:6:end,veh,:), time_horizon+1, mc_samples),[],2)-x_mean;
    yneg = min(reshape(x_mean_holder(3:6:end,veh,:), time_horizon+1, mc_samples),[],2)-y_mean;
    ypos = max(reshape(x_mean_holder(3:6:end,veh,:), time_horizon+1, mc_samples),[],2)-y_mean;
    errorbar(x_mean,y_mean,yneg,ypos,xneg,xpos,'--', 'Color', colors(veh,:))
    patch('Faces',F_xz,'Vertices', Polyhedron(target_sets(veh).A([1;2;3;7;8;9],1:3), target_sets(veh).b([1;2;3;7;8;9])).V(:,[1,3]),...
        'FaceColor', colors(veh, :), ...
        'FaceAlpha', 0.1); 
    plot(x_0_not_random(1,veh), x_0_not_random(3,veh), 'Color', 'k', 'Marker', plot_symbols(veh), 'MarkerFaceColor', 'k', 'LineStyle','none');
end
xlabel('x (in meters)')
ylabel('z (in meters)')
hold off