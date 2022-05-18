%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make some plots
%%%%%%%%%%%%%%%%%%%%%%%%%%


colors = [224,   0,   0; % red
           30, 144, 255; % dark blue
           0, 170,  85   % green
           ] ./ 255;
       
plot_symbols = ['o', 'd'];
F_xy = [1,2,3,4];

dcm = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];

fig = figure();
fig.Units    = 'inches';
fig.Position = [1,1,19,8];

subplot(7,2,[1 2]);
hold on 
plot(nan, nan, 'Color', colors(2, :), 'Marker', plot_symbols(1));
plot(nan, nan, 'Color', colors(3, :), 'Marker', plot_symbols(2));
quiver(x_0_deputy(1),x_0_deputy(2),0, -1, 0, 'Color', 'k');

patch('Faces',F_xy,'Vertices', Polyhedron([-eye(2); eye(2)], ones(4,1)).V,...
    'FaceColor',  'white', ...
    'EdgeColor', 'black', ...
    'FaceAlpha', 0); 


p = nsidedpoly(1000, 'Center', [100 0], 'Radius', r);
plot(p, 'EdgeColor', colors(1,:), 'FaceColor', colors(1,:), 'LineStyle','--', 'FaceAlpha', 0.1)

plot(x_0_deputy(1), nan, 'Color', 'k', 'Marker', 's', 'MarkerFaceColor', 'k', 'LineStyle','none');
plots=get(gca, 'Children');

legend([plots(6), plots(5), plots(4), plots(3), plots(2), plots(1)], ...
     {'Proposed Method', 'Particle Control', 'Sensor Direction', 'Target Set', 'Collision Avoidance Region', 'Initial Location' },...
    'Orientation','horizontal', ...
    'Location', 'south', ...
    'NumColumns', 3, ...
    'interpreter', 'latex');


axis([0 0.1 0 0.1]);
axis off
hold off


subplot(7,2,[3:2:14]);
hold on

p = nsidedpoly(1000, 'Center', [0 0], 'Radius', r);
plot(p, 'EdgeColor', colors(1,:), 'FaceColor', colors(1,:), 'LineStyle','--', 'FaceAlpha', 0.1)

patch('Faces',F_xy,'Vertices', Polyhedron(target_set_2.A([1;2;7;8],1:2), target_set_2.b([1;2;7;8])).V,...
    'FaceColor', 'white', ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_4.A([1;2;7;8],1:2), target_set_4.b([1;2;7;8])).V,...
    'FaceColor', 'white', ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_6.A([1;2;7;8],1:2), target_set_6.b([1;2;7;8])).V,...
    'FaceColor', 'white', ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_8.A([1;2;7;8],1:2), target_set_8.b([1;2;7;8])).V,...
    'FaceColor', 'white', ...
    'FaceAlpha', 0.1); 

for i = 2:2:time_horizon
   index = 6 * (i-1) + [1:3];
   mu_i = mean_X(index);
   theta = -1/2*pi - i/8*pi;
   vec_upper  = dcm(theta+theta_range)*[0;-1];
   vec_lower  = dcm(theta-theta_range)*[0;-1];
   p=quiver(mu_i(1), mu_i(2), vec_lower(1), vec_lower(2), 0, 'Color', 'k');
   p.MaxHeadSize=0;
   p=quiver(mu_i(1), mu_i(2), vec_upper(1), vec_upper(2), 0, 'Color', 'k');
   p.MaxHeadSize=0;
end

plot([x_0_deputy(1); mean_X(1:6:end)], [x_0_deputy(2); mean_X(2:6:end)], 'Color', colors(2, :), 'Marker', plot_symbols(1));
plot(x_0_deputy(1), x_0_deputy(2), 'Color', 'k', 'Marker', 's', 'MarkerFaceColor', 'k', 'LineStyle','none');

for i = 1:time_horizon
   index = 6 * (i-1) + [1:3];
   mu_i = mean_X(index);
   vec  = dcm(mu_i(3))*[0;-2];
   p=quiver(mu_i(1), mu_i(2), vec(1), vec(2), 0, 'Color', colors(2, :));
   p.MaxHeadSize=5;
end

p=quiver(x_0_deputy(1),x_0_deputy(2),0, -2, 0, 'Color', 'k');
p.MaxHeadSize=5;


xlabel('x (in meters)')
ylabel('y (in meters)')
 
axis([-14 12 -12 4])
axis equal
hold off

subplot(7,2,[4:2:14]);
hold on


p = nsidedpoly(8, 'Center', [0 0], 'Radius', 8.659);
plot(p, 'EdgeColor', colors(1,:), 'FaceColor', colors(1,:), 'LineStyle','--', 'FaceAlpha', 0.1)

patch('Faces',F_xy,'Vertices', Polyhedron(target_set_2.A([1;2;7;8],1:2), target_set_2.b([1;2;7;8])).V,...
    'FaceColor', 'white', ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_4.A([1;2;7;8],1:2), target_set_4.b([1;2;7;8])).V,...
    'FaceColor', 'white', ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_6.A([1;2;7;8],1:2), target_set_6.b([1;2;7;8])).V,...
    'FaceColor', 'white', ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_8.A([1;2;7;8],1:2), target_set_8.b([1;2;7;8])).V,...
    'FaceColor', 'white', ...
    'FaceAlpha', 0.1); 

for i = 2:2:time_horizon
   index = 6 * (i-1) + [1:3];
   mu_i = mean_X_pc(index);
   theta = -1/2*pi - i/8*pi;
   vec_upper  = dcm(theta+theta_range)*[0;-1];
   vec_lower  = dcm(theta-theta_range)*[0;-1];
   p=quiver(mu_i(1), mu_i(2), vec_lower(1), vec_lower(2), 0, 'Color', 'k');
   p.MaxHeadSize=0;
   p=quiver(mu_i(1), mu_i(2), vec_upper(1), vec_upper(2), 0, 'Color', 'k');
   p.MaxHeadSize=0;
end

plot([x_0_deputy(1); mean_X_pc(1:6:end)], [x_0_deputy(2); mean_X_pc(2:6:end)], 'Color', colors(3, :), 'Marker', plot_symbols(2));

for i = 1:time_horizon
   index = 6 * (i-1) + [1:3];
   mu_i = mean_X_pc(index);
   vec  = dcm(mu_i(3))*[0;-2];
   p=quiver(mu_i(1), mu_i(2), vec(1), vec(2), 0, 'Color', colors(3, :));
   p.MaxHeadSize=5;
end

plot(x_0_deputy(1), x_0_deputy(2), 'Color', 'k', 'Marker', 's', 'MarkerFaceColor', 'k', 'LineStyle','none');
p=quiver(x_0_deputy(1),x_0_deputy(2),0, -2, 0, 'Color', 'k');
p.MaxHeadSize=5;



xlabel('x (in meters)')
ylabel('y (in meters)')
 
axis([-14 12 -12 4])
axis equal
hold off
