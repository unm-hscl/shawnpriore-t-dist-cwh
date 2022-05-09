%% set up satellites
% universal
sc_XR=[0 0 0 1 0 0 0 ]';     
t0=[2020 01 15 00 00 00];
drag_coefficient = 0.2;


% Set chief sat
x_0_chief = [orbital_radius; 0; 0; 0; sqrt(gravitational_body/ orbital_radius); 0];
chief_craft_mass      = 3000;     % kg
chief_area_exposed    = 5^2;      % m^2
chief_reflectivity    = 0.25;     % unitless
sc_chief = InertialSpacecraftModel(x_0_chief, sc_XR, t0, chief_craft_mass, chief_area_exposed, chief_reflectivity, drag_coefficient);


% docking sat params
other_craft_mass      = 100;       % kg
other_area_exposed    = 0.5^2;     % m^2
other_reflectivity    = 0.3;       % unitless

% docking sat init conditions
x_0_a_inertial = x_0_chief + [x_0_a(2); -x_0_a(1); zeros(4,1)];
x_0_b_inertial = x_0_chief + [x_0_b(2); -x_0_b(1); zeros(4,1)];
x_0_c_inertial = x_0_chief + [x_0_c(2); -x_0_c(1); zeros(4,1)];
x_0_d_inertial = x_0_chief + [x_0_d(2); -x_0_d(1); zeros(4,1)];
x_0_e_inertial = x_0_chief + [x_0_e(2); -x_0_e(1); zeros(4,1)];

sc_a = InertialSpacecraftModel(x_0_a_inertial, sc_XR, t0, other_craft_mass, other_area_exposed, other_reflectivity, drag_coefficient);
sc_b = InertialSpacecraftModel(x_0_b_inertial, sc_XR, t0, other_craft_mass, other_area_exposed, other_reflectivity, drag_coefficient);
sc_c = InertialSpacecraftModel(x_0_c_inertial, sc_XR, t0, other_craft_mass, other_area_exposed, other_reflectivity, drag_coefficient);
sc_d = InertialSpacecraftModel(x_0_d_inertial, sc_XR, t0, other_craft_mass, other_area_exposed, other_reflectivity, drag_coefficient);
sc_e = InertialSpacecraftModel(x_0_e_inertial, sc_XR, t0, other_craft_mass, other_area_exposed, other_reflectivity, drag_coefficient);



% adjust controller 
dcm = zeros(time_horizon*3);

for i = 0:(time_horizon-1)
    index = i*3 + [1:3];
    t = i*sampling_period;
    dcm(index, index) = [cos(pi/2-orbit_ang_vel*t) sin(pi/2-orbit_ang_vel*t) 0; 
                        -sin(pi/2-orbit_ang_vel*t) cos(pi/2-orbit_ang_vel*t) 0; 
                        0 0 1];
end

U_a_rot = dcm*U_a;
U_b_rot = dcm*U_b;
U_c_rot = dcm*U_c;
U_d_rot = dcm*U_d;
U_e_rot = dcm*U_e;


% propigate motion
for i = 1:time_horizon
    index = (i-1)*3 + [1:3];
    
    sc_chief.PropTransEOM(sampling_period);
    
    sc_a.XT_traj(end,4:6) = sc_a.XT_traj(end,4:6) + U_a_rot(index)' * sc_a.v_scale;
    sc_a.PropTransEOM(sampling_period);
    
    sc_b.XT_traj(end,4:6) = sc_b.XT_traj(end,4:6) + U_b_rot(index)' * sc_b.v_scale;
    sc_b.PropTransEOM(sampling_period);
    
    sc_c.XT_traj(end,4:6) = sc_c.XT_traj(end,4:6) + U_c_rot(index)' * sc_c.v_scale;
    sc_c.PropTransEOM(sampling_period);
    
    sc_d.XT_traj(end,4:6) = sc_d.XT_traj(end,4:6) + U_d_rot(index)' * sc_d.v_scale;
    sc_d.PropTransEOM(sampling_period);
    
    sc_e.XT_traj(end,4:6) = sc_e.XT_traj(end,4:6) + U_e_rot(index)' * sc_e.v_scale;
    sc_e.PropTransEOM(sampling_period);

end


%% Create relative motion object
rel_a = RelativeMotion(sc_chief, sc_a);
rel_b = RelativeMotion(sc_chief, sc_b);
rel_c = RelativeMotion(sc_chief, sc_c);
rel_d = RelativeMotion(sc_chief, sc_d);
rel_e = RelativeMotion(sc_chief, sc_e);

rel_a_XT_Rel_traj_cwh = zeros(size(rel_a.XT_Rel_traj));
rel_b_XT_Rel_traj_cwh = zeros(size(rel_b.XT_Rel_traj));
rel_c_XT_Rel_traj_cwh = zeros(size(rel_c.XT_Rel_traj));
rel_d_XT_Rel_traj_cwh = zeros(size(rel_d.XT_Rel_traj));
rel_e_XT_Rel_traj_cwh = zeros(size(rel_e.XT_Rel_traj));

% rotate back to cwh for comparison
for i = 1:size(rel_a.T_traj,1)
    t = rel_a.T_traj(i)-rel_a.T_traj(1);
    dcm = [cos(pi/2-orbit_ang_vel*t) -sin(pi/2-orbit_ang_vel*t) 0; 
           sin(pi/2-orbit_ang_vel*t) cos(pi/2-orbit_ang_vel*t) 0; 
           0 0 1];
    rel_a_XT_Rel_traj_cwh(i,1:3) = (dcm * rel_a.XT_Rel_traj(i,1:3)')';
    rel_b_XT_Rel_traj_cwh(i,1:3) = (dcm * rel_b.XT_Rel_traj(i,1:3)')';
    rel_c_XT_Rel_traj_cwh(i,1:3) = (dcm * rel_c.XT_Rel_traj(i,1:3)')';
    rel_d_XT_Rel_traj_cwh(i,1:3) = (dcm * rel_d.XT_Rel_traj(i,1:3)')';
    rel_e_XT_Rel_traj_cwh(i,1:3) = (dcm * rel_e.XT_Rel_traj(i,1:3)')';
end

% make plots
figure()
hold on
pa = plot3(rel_a_XT_Rel_traj_cwh(:,1), rel_a_XT_Rel_traj_cwh(:,2), rel_a_XT_Rel_traj_cwh(:,3), 'Color', red, 'Marker', 'h', 'LineWidth', 1);
pb = plot3(rel_b_XT_Rel_traj_cwh(:,1), rel_b_XT_Rel_traj_cwh(:,2), rel_b_XT_Rel_traj_cwh(:,3), 'Color', blue, 'Marker', 'p', 'LineWidth', 1);
pc = plot3(rel_c_XT_Rel_traj_cwh(:,1), rel_c_XT_Rel_traj_cwh(:,2), rel_c_XT_Rel_traj_cwh(:,3), 'Color', green, 'Marker', '^', 'LineWidth', 1);
pd = plot3(rel_d_XT_Rel_traj_cwh(:,1), rel_d_XT_Rel_traj_cwh(:,2), rel_d_XT_Rel_traj_cwh(:,3), 'Color', purple, 'Marker', 'o', 'LineWidth', 1);
pe = plot3(rel_e_XT_Rel_traj_cwh(:,1), rel_e_XT_Rel_traj_cwh(:,2), rel_e_XT_Rel_traj_cwh(:,3), 'Color', grey, 'Marker', 's', 'LineWidth', 1);
p1 = patch('Faces',F,'Vertices', Polyhedron(target_set_a.A([1;2;3;7;8;9],1:3), target_set_a.b([1;2;3;7;8;9])).V,...
    'FaceColor', red, ...
    'FaceAlpha',0.1); 
patch('Faces',F,'Vertices', Polyhedron(target_set_b.A([1;2;3;7;8;9],1:3), target_set_b.b([1;2;3;7;8;9])).V,...
    'FaceColor', blue, ...
    'FaceAlpha',0.1); 
patch('Faces',F,'Vertices', Polyhedron(target_set_c.A([1;2;3;7;8;9],1:3), target_set_c.b([1;2;3;7;8;9])).V,...
    'FaceColor', green, ...
    'FaceAlpha',0.1); 
patch('Faces',F,'Vertices', Polyhedron(target_set_d.A([1;2;3;7;8;9],1:3), target_set_d.b([1;2;3;7;8;9])).V,...
    'FaceColor', purple, ...
    'FaceAlpha',0.1); 
patch('Faces',F,'Vertices', Polyhedron(target_set_e.A([1;2;3;7;8;9],1:3), target_set_e.b([1;2;3;7;8;9])).V,...
    'FaceColor', grey, ...
    'FaceAlpha',0.1); 
p2 = plot3(x_0_a(1), x_0_a(2), x_0_a(3), 'Color', 'k', 'Marker', 'h', 'MarkerFaceColor', 'k', 'LineStyle','none');
plot3(x_0_b(1), x_0_b(2), x_0_b(3), 'Color', 'k', 'Marker', 'p', 'MarkerFaceColor', 'k');
plot3(x_0_c(1), x_0_c(2), x_0_c(3), 'Color', 'k', 'Marker', '^', 'MarkerFaceColor', 'k');
plot3(x_0_d(1), x_0_d(2), x_0_d(3), 'Color', 'k', 'Marker', 'o', 'MarkerFaceColor', 'k');
plot3(x_0_e(1), x_0_e(2), x_0_e(3), 'Color', 'k', 'Marker', 's', 'MarkerFaceColor', 'k');

xlabel('x (in meters)')
ylabel('y (in meters)')
zlabel('z (in meters)')

l = legend([pa,pb,pc,pd,pe,p1,p2], {'A', 'B', 'C', 'D', 'E', 'Target Set', 'Initial Location' },...
    'Orientation','horizontal', ...
    'Location', 'northoutside', ...
    'NumColumns', 4, ...
    'interpreter', 'latex');

set(gca, 'OuterPosition', [0.025,(1-0.95/2)/2,0.95,0.95/2]);
hold off
drawnow()
view(3)


