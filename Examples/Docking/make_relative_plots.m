%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make some plots
%%%%%%%%%%%%%%%%%%%%%%%%%%

F_xy = [1, 4, 6, 8];
F_xz = [3, 4, 6, 5];

colors = [224,   0,   0; % red
           30, 144  255; % dark blue
            0, 170,  85; % green
          118,   0, 168; % purple
           46,  52,  59; % grey
          236, 176,  31;  % yellow
           76, 189, 237; % light blue
          161,  19,  46  % dark red
           ] ./ 255;
       
plot_symbols = ['h', 'p', '^', 'o', 'v', '>', 'd', 's'];

fig = figure();
fig.Units    = 'inches';
fig.Position = [1,1,19,8];

subplot(7,2,[1 2]);
hold on 
for i = 1:vehicles
    plot(nan, nan, 'Color', colors(i, :), 'Marker', plot_symbols(i));
end
patch('Faces',F_xy,'Vertices', Polyhedron(target_sets(1).A([1;2;3;7;8;9],1:3), target_sets(1).b([1;2;3;7;8;9])).V,...
    'FaceColor',  'white', ...
    'EdgeColor', 'black', ...
    'FaceAlpha', 0); 
plot(x_0(1,1), nan, 'Color', 'k', 'Marker', plot_symbols(1), 'MarkerFaceColor', 'k', 'LineStyle','none');
plots=get(gca, 'Children');

legend([plots(9), plots(8), plots(7), plots(6), plots(5), plots(4), plots(3), plots(2), plots(1)], ...
     {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'Target Set', 'Initial Location' },...
    'Orientation','horizontal', ...
    'Location', 'south', ...
    'NumColumns', 9, ...
    'interpreter', 'latex');


axis([0 1 0 1]);
axis off
hold off


subplot(7,2,[3:2:14]);
hold on
for i = 1:vehicles
    plot([x_0(1,i); mean_X(1:6:end,i)], [x_0(2,i); mean_X(2:6:end,i)], 'Color', colors(i, :), 'Marker', plot_symbols(i));
    patch('Faces',F_xy,'Vertices', Polyhedron(target_sets(i).A([1;2;3;7;8;9],1:3), target_sets(i).b([1;2;3;7;8;9])).V,...
        'FaceColor', colors(i, :), ...
        'FaceAlpha', 0.1); 
    plot(x_0(1,i), x_0(2,i), 'Color', 'k', 'Marker', plot_symbols(i), 'MarkerFaceColor', 'k', 'LineStyle','none');
end
xlabel('x (in meters)')
ylabel('y (in meters)')
hold off

subplot(7,2,[4:2:14]);
hold on
for i = 1:vehicles
    
    plot([x_0(1,i); mean_X(1:6:end,i)], [x_0(3,i); mean_X(3:6:end,i)], 'Color', colors(i, :), 'Marker', plot_symbols(i));
    patch('Faces',F_xz,'Vertices', Polyhedron(target_sets(i).A([1;2;3;7;8;9],1:3), target_sets(i).b([1;2;3;7;8;9])).V(:,[1,3]),...
        'FaceColor', colors(i, :), ...
        'FaceAlpha', 0.1); 
    plot(x_0(1,i), x_0(3,i), 'Color', 'k', 'Marker', plot_symbols(i), 'MarkerFaceColor', 'k', 'LineStyle','none');
end
xlabel('x (in meters)')
ylabel('z (in meters)')
hold off