%% inputs:

% clear
% clc

% model parameters:
stiffness_g_t = 100; % ground stiffness in tension (must be > 0)
stiffness_g_c = 10; % ground stiffness in compression (must be > 0)
stiffness_f = 20; % fiber stiffness (must be > 0)
straightness_f = 0.83; % fiber straightness (must be > 0 and <= 1)
length_sample_axial = 200;
length_sample_trans = 100;
num_seeds = 200;
disp_x1_prescribed = 15; % magnitude of x1 displacement at the right boundary for each step
steps = 60;


%% initialization -- generate the network:

% preset general function settings:
sympref('HeavisideAtOrigin', 1);
solver_options = optimset('Display','off');

% generate network from random seed points:
x1_range = [0 length_sample_axial+1]; % x1, axial dimension;
x2_range = [0 length_sample_trans+1]; % x2, transverse dimension;
x1 = randi(x1_range, 1, num_seeds);
x2 = randi(x2_range, 1, num_seeds);
[~,~,vx1,vx2] = VoronoiLimitRectSquare(x1, x2, [x1_range x2_range], 1, 1); % Author: Preetham Manjunatha (https://github.com/preethamam)
vx1 = vx1 - 1; vx2 = vx2 - 1; % VoronoiLimitRectSquare() only gives vertex positions >= 1, but we want >= 0, so we +1 to range max then -1 to all vertex positions
hold on

% record number of edges (~fibers) and vertices (~crosslinks):
vx1_noSelfEdges = [];
vx2_noSelfEdges = [];
for i_edge = 1:size(vx1,2)
    if (vx1(1,i_edge) ~= vx1(2,i_edge) || vx2(1,i_edge) ~= vx2(2,i_edge))
        vx1_noSelfEdges(:,end+1) = vx1(:,i_edge);
        vx2_noSelfEdges(:,end+1) = vx2(:,i_edge);
    end
end
num_edges = size(vx1_noSelfEdges,2);
verts_all = [vx1_noSelfEdges(1,:) vx1_noSelfEdges(2,:); vx2_noSelfEdges(1,:) vx2_noSelfEdges(2,:)]';
[verts_unique,~,i_unique] = unique(verts_all, 'rows'); % i.e. verts_all = verts_unique(i_unique,:)
num_verts = size(verts_unique, 1);

% construct the reference connectivity matrix:
con_verts = zeros(num_verts); % reference connectivity matrix; symmetrical
connections = [i_unique(1:num_edges) i_unique(num_edges+1:end)]; % pairs of indices in each row represent connected vertices
for i = 1:num_edges
    con_verts(connections(i,1),connections(i,2)) = 1;
    con_verts(connections(i,2),connections(i,1)) = 1; % because matrix should be symmetrical
end
con_verts_old = con_verts;

% store vertex positions:
pos_verts_ref = [verts_unique(:,1); verts_unique(:,2)]; % reference position of each vertex, stored as [x1.1,x1.2,...x1.n,x2.1,x2.2,...x2.n]'

% group vertices based on position:
i_verts_group_L = find(pos_verts_ref(1:num_verts) == 0); % vertices at left boundary (fixed)
i_verts_group_R = find(pos_verts_ref(1:num_verts) == x1_range(2)-1); % vertices at right boundary (prescribed displacement)

% remove middle vertices that have only a single neighbor:
i_verts_group_singles = find(sum(con_verts) <= 1); % all vertices with a single neighbor
[Lia,Locb] = ismember([i_verts_group_L' i_verts_group_R'],i_verts_group_singles);
i_verts_group_singles(Locb(Lia)) = []; % exclude L and R boundary vertices from this group
con_verts(i_verts_group_singles,:) = []; con_verts(:,i_verts_group_singles) = []; % remove middle singles from connectivity matrix
pos_verts_ref([i_verts_group_singles' i_verts_group_singles'+num_verts],:) = []; % remove middle singles from positions array
num_verts = size(con_verts,1); % account for the fact that number of vertices has changed (decreased)

% re-group vertices based on position (because con_verts_ref and pos_verts_ref have changed size):
i_verts_group_L = find(pos_verts_ref(1:num_verts) == 0); % vertices at left boundary (fixed)
i_verts_group_R = find(pos_verts_ref(1:num_verts) == x1_range(2)-1); % vertices at right boundary (prescribed displacement)
i_verts_group_M = find(pos_verts_ref(1:num_verts) > 0 & pos_verts_ref(1:num_verts) < x1_range(2)-1); % vertices in middle (to be equilibrated)

% randomly offset reference vertex positions:
pos_verts_ref(num_verts+1:end) = pos_verts_ref(num_verts+1:end) + randn(num_verts,1)/10; % x2 positions, for all vertices
pos_verts_ref(i_verts_group_M) = pos_verts_ref(i_verts_group_M) + randn(length(i_verts_group_M),1)/10; % x1 positions, for middle vertices only

% scatter(pos_verts_ref(i_verts_group_L)+1, pos_verts_ref(i_verts_group_L+num_verts)+1, 'g')
% scatter(pos_verts_ref(i_verts_group_M)+1, pos_verts_ref(i_verts_group_M+num_verts)+1, 'k')
% scatter(pos_verts_ref(i_verts_group_R)+1, pos_verts_ref(i_verts_group_R+num_verts)+1, 'r')


%% initialize other arrays:
pos_verts_cur = pos_verts_ref; % current position of each vertex
pos_verts_allSteps = zeros(steps, num_verts, 2); % contains equilibrium positions for all steps, stored as [step, vertex, coordinate (x1:x2)]
F1_verts_R_allSteps = zeros(steps, length(i_verts_group_R)); % x1 force at each R boundary vertex for all steps, stored as [step, vertex]
F_net_residual_allSteps = zeros(steps, length(i_verts_group_M)); % net force remaining after equlibration at each middle vertex for all steps, stored as [step, vertex]


%% computation -- iteratively solve the model:

% iterate through displacement steps:
for i_step = 1:steps
    
    % apply the prescribed x1 displacement linearly across the vertices (i.e. max at R boundary to 0 at L boundary):
    if i_step > 1
        disp_linear = zeros(num_verts*2, 1);
        disp_linear(1:num_verts) = disp_x1_prescribed / max(pos_verts_cur(i_verts_group_R)) .* pos_verts_cur(1:num_verts);
        pos_verts_cur = pos_verts_cur + disp_linear; % naively update x1 positions (these will be used as initial guesses for the solver below)
    end

    % loop repeatedly through middle vertices until all forces approach 0:
    F_net_residual = ones(length(i_verts_group_M),1);
    F_net_residual_prev = zeros(length(i_verts_group_M),1);
    num_loopsCompleted = 0;
    while sqrt(mean((F_net_residual-F_net_residual_prev).^2)) > 0.1 % use RMS error between loops
        % store previous residual net forces:
        F_net_residual_prev = F_net_residual;
        % loop once through middle vertices to equilibrate their positions sequentially:
        for i_vert = 1:length(i_verts_group_M)
            % vertex of interest:
            X = [pos_verts_ref(i_verts_group_M(i_vert)), pos_verts_ref(i_verts_group_M(i_vert) + num_verts)]; % vertex's reference position
            x_guess = [pos_verts_cur(i_verts_group_M(i_vert)), pos_verts_cur(i_verts_group_M(i_vert) + num_verts)]; % starting guess for vertex's current position
            % neighbor vertices:
            i_connections = find(con_verts(i_verts_group_M(i_vert),:)); % 1D array containing indices of all the vertex's connections
            V = [pos_verts_ref(i_connections), pos_verts_ref(i_connections + num_verts)]; % neighbors' reference positions, stored as [vertex, coordinate (x1:x2)]
            v = [pos_verts_cur(i_connections), pos_verts_cur(i_connections + num_verts)]; % neighbors' current positions, stored as [vertex, coordinate (x1:x2)]
            % system of equations:
            F_net_residual(i_vert) = norm(forceFunctions(X,x_guess,V,v,stiffness_f,straightness_f,stiffness_g_t,stiffness_g_c)); % store residual net force (before equilibration)
            system = @(x)forceFunctions(X,x,V,v,stiffness_f,straightness_f,stiffness_g_t,stiffness_g_c); % functions to solve (see "helper functions" below)
            x_eq = fsolve(system, x_guess, solver_options); % solve system to find vertex's current equilibrium position
            % update vertex position:
            pos_verts_cur(i_verts_group_M(i_vert)) = x_eq(1); % update x1
            pos_verts_cur(i_verts_group_M(i_vert) + num_verts) = x_eq(2); % update x2
        end
        % display progress:
        num_loopsCompleted = num_loopsCompleted + 1;
        disp(['loop complete: step, ' num2str(i_step) ' | loop, ' num2str(num_loopsCompleted) ' | mean residual force, ' num2str(mean(F_net_residual))]);
    end

    % store equilibrated vertex positions for this step:
    pos_verts_allSteps(i_step,:,1) = pos_verts_cur(1:num_verts); % x1
    pos_verts_allSteps(i_step,:,2) = pos_verts_cur(num_verts+1:end); % x2
    
    % store x1 forces at R boundary vertices for this step:
    for i_vert = 1:length(i_verts_group_R)
        % vertex of interest:
        X = [pos_verts_ref(i_verts_group_R(i_vert)), pos_verts_ref(i_verts_group_R(i_vert) + num_verts)]; % vertex's reference position
        x = [pos_verts_cur(i_verts_group_R(i_vert)), pos_verts_cur(i_verts_group_R(i_vert) + num_verts)]; % vertex's current position
        % neighbor vertices:
        i_connections = find(con_verts(i_verts_group_R(i_vert),:)); % 1D array containing indices of all the vertex's connections
        V = [pos_verts_ref(i_connections), pos_verts_ref(i_connections + num_verts)]; % neighbors' reference positions, stored as [vertex, coordinate (x1:x2)]
        v = [pos_verts_cur(i_connections), pos_verts_cur(i_connections + num_verts)]; % neighbors' current positions, stored as [vertex, coordinate (x1:x2)]
        % calculate force:
        F_verts_R = forceFunctions(X,x,V,v,stiffness_f,straightness_f,stiffness_g_t,stiffness_g_c);
        F1_verts_R_allSteps(i_step, i_vert) = F_verts_R(1); % store only x1 component of force
    end

    % store remaining net residual forces at middle vertices for this step:
    F_net_residual_allSteps(i_step,:) = F_net_residual;
    
end

% visualize:
figure
plot(1:size(F1_verts_R_allSteps,1), -F1_verts_R_allSteps)
title('R boundary forces')

figure
plot(1:size(F_net_residual_allSteps,1), F_net_residual_allSteps)
title('M residual forces')
ylim([0,max(max(abs(F1_verts_R_allSteps)))])

figure
for i = 1:size(pos_verts_allSteps, 1)
    scatter(squeeze(pos_verts_allSteps(i,:,1)), squeeze(pos_verts_allSteps(i,:,2)), 15, ...
        'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [i/size(pos_verts_allSteps, 1) i/size(pos_verts_allSteps, 1) 0.5], ...
        'LineWidth', 0.3); hold on
end
title('displacements')
% draw reference outline:
x_max = max(pos_verts_ref(1:num_verts));
y_max = max(pos_verts_ref(num_verts+1:end));
y_min = min(pos_verts_ref(num_verts+1:end));
line([x_max x_max], [y_min y_max], 'Color', [1 0 0], 'LineStyle','--', 'LineWidth', 1);
% draw current outline:
x_max = max(pos_verts_cur(1:num_verts));
line([x_max x_max], [y_min y_max], 'Color', [0 1 0], 'LineStyle','--', 'LineWidth', 1);


%% data processing:

% calculate axial and transverse stretchs:
stretch_axial = max(pos_verts_allSteps(:,:,1),[],2) / max(pos_verts_allSteps(1,:,1));
pos_verts_x1center_x2minMax = zeros(steps, 2);
for i_step = 1:steps
    pos_verts_x1center_x2minMax(i_step,1) = min(pos_verts_allSteps(i_step, pos_verts_allSteps(i_step,:,1) > max(pos_verts_allSteps(i_step,:,1))/2.2 ...
        & pos_verts_allSteps(i_step,:,1) < max(pos_verts_allSteps(i_step,:,1))/1.8, 2));
    pos_verts_x1center_x2minMax(i_step,2) = max(pos_verts_allSteps(i_step, pos_verts_allSteps(i_step,:,1) > max(pos_verts_allSteps(i_step,:,1))/2.2 ...
        & pos_verts_allSteps(i_step,:,1) < max(pos_verts_allSteps(i_step,:,1))/1.8, 2));
end
stretch_trans = (pos_verts_x1center_x2minMax(:,2) - pos_verts_x1center_x2minMax(:,1)) / (pos_verts_x1center_x2minMax(1,2) - pos_verts_x1center_x2minMax(1,1));

% plot axial stretch-transverse stretch and axial stretch-force relations:
figure
subplot(1,2,1)
plot(stretch_axial, stretch_trans)
xlabel('axial stretch')
ylabel('transverse stretch')
ylim([0 1])
xlim([1 4.5])
subplot(1,2,2)
plot(stretch_axial, -sum(F1_verts_R_allSteps,2))
ylim([0 20000])
xlim([1 4.5])
xlabel('axial stretch')
ylabel('mean F1 force')


%% visual checks:

i_step = 20;

% show edge tensions:
subplot(1,2,1)
scatter(squeeze(pos_verts_allSteps(i_step,:,1)), squeeze(pos_verts_allSteps(i_step,:,2)), 1, ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'LineWidth', 0.1); hold on
for i_vertA = 1:size(con_verts,1)
    for i_vertB = 1:size(con_verts,2)
        if con_verts(i_vertA,i_vertB) == 1
            color_fiber = [0.7 0.7 0.7];
            length_fiber_ref = sqrt((pos_verts_allSteps(1,i_vertA,1)-pos_verts_allSteps(1,i_vertB,1)).^2 ...
                + (pos_verts_allSteps(1,i_vertA,2)-pos_verts_allSteps(1,i_vertB,2)).^2);
            length_fiber_cur = sqrt((pos_verts_allSteps(i_step,i_vertA,1)-pos_verts_allSteps(i_step,i_vertB,1)).^2 ...
                + (pos_verts_allSteps(i_step,i_vertA,2)-pos_verts_allSteps(i_step,i_vertB,2)).^2);
            if length_fiber_cur/length_fiber_ref > 1/straightness_f % fully loaded in tension
                color_fiber = [1 0 0];
            elseif length_fiber_cur/length_fiber_ref > 1 % only ground loaded in tension
                color_fiber = [0.8 0.7 0];
            elseif length_fiber_cur/length_fiber_ref < 1 % only ground loaded in compression
                color_fiber = [0 0 1];
            end
            plot([pos_verts_allSteps(i_step,i_vertA,1) pos_verts_allSteps(i_step,i_vertB,1)], ...
                [pos_verts_allSteps(i_step,i_vertA,2) pos_verts_allSteps(i_step,i_vertB,2)], ...
                'Color', color_fiber, 'LineWidth', 2);
        end
    end
end
ylim([0 100]);

% show edge reorientations:
subplot(1,2,2)
scatter(squeeze(pos_verts_allSteps(i_step,:,1)), squeeze(pos_verts_allSteps(i_step,:,2)), 1, ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'LineWidth', 0.1); hold on
for i_vertA = 1:size(con_verts,1)
    for i_vertB = 1:size(con_verts,2)
        if con_verts(i_vertA,i_vertB) == 1
            angle_fiber_ref = atan2(abs(pos_verts_allSteps(1,i_vertA,2)-pos_verts_allSteps(1,i_vertB,2)), ...
                abs(pos_verts_allSteps(1,i_vertA,1)-pos_verts_allSteps(1,i_vertB,1)));
            angle_fiber_cur = atan2(abs(pos_verts_allSteps(i_step,i_vertA,2)-pos_verts_allSteps(i_step,i_vertB,2)), ...
                abs(pos_verts_allSteps(i_step,i_vertA,1)-pos_verts_allSteps(i_step,i_vertB,1)));
            color_fiber = [1-angle_fiber_cur/(pi/2) 0 angle_fiber_cur/(pi/2)];
            plot([pos_verts_allSteps(i_step,i_vertA,1) pos_verts_allSteps(i_step,i_vertB,1)], ...
                [pos_verts_allSteps(i_step,i_vertA,2) pos_verts_allSteps(i_step,i_vertB,2)], ...
                'Color', color_fiber, 'LineWidth', 2);
        end
    end
end
ylim([0 100]);




%% helper functions:

function F = forceFunctions(X,x,V,v,k_f,s_f,k_g_t,k_g_c)

% X(1:2) is the known reference position vector of the vertex of interest (VOI)
% x(1:2) is the unknown current position vector of the VOI ==> this is what we're solving for!
% V(:,1:2) is an array containing the reference position vectors of all the vertices our VOI is connected to
% v(:,1:2) is an array containing the current position vectors of all the vertices our VOI is connected to

% initialize force components:
F(1) = 0;
F(2) = 0;

% add forces from neighboring vertices:
for n = 1:size(V,1)
    l = sqrt((v(n,1)-x(1))^2 + (v(n,2)-x(2))^2); % current edge length
    L_g = sqrt((V(n,1)-X(1))^2 + (V(n,2)-X(2))^2); % reference edge length, ground
    L_f = L_g / s_f; % reference edge length, fiber
    if l/L_f > 1
        % include fibers and ground in force contribution:
        if 1/L_g >= 1
            % use tensile (t) stiffness for ground:
            F(1) = F(1) + (k_f*(l-L_f)^2 + k_g_t*(l-L_g)) * (v(n,1)-x(1))/l; % linear: k_f * (l-L_f)
            F(2) = F(2) + (k_f*(l-L_f)^2 + k_g_t*(l-L_g)) * (v(n,2)-x(2))/l;
        else
            % use compressive (c) stiffness for ground:
            F(1) = F(1) + (k_f*(l-L_f)^2 + k_g_c*(l-L_g)) * (v(n,1)-x(1))/l;
            F(2) = F(2) + (k_f*(l-L_f)^2 + k_g_c*(l-L_g)) * (v(n,2)-x(2))/l;
        end
    else
        % include only ground in force contribution:
        if 1/L_g >= 1
            % use tensile (t) stiffness for ground:
            F(1) = F(1) + k_g_t*(l-L_g) * (v(n,1)-x(1))/l;
            F(2) = F(2) + k_g_t*(l-L_g) * (v(n,2)-x(2))/l;
        else
            % use compressive (c) stiffness for ground:
            F(1) = F(1) + k_g_c*(l-L_g) * (v(n,1)-x(1))/l;
            F(2) = F(2) + k_g_c*(l-L_g) * (v(n,2)-x(2))/l;
        end
    end
end

end



%% trashed code:

% X = [1 1]; % vertex reference position
% x = [1 1]; % starting guess for vertex current position
% V = [0.5 3; 2 2; 4 0.5; 0.5 0]; % neighbor reference positions
% v = [0.5 3; 2.5 2; 5 0.5; 0.5 0]; % neighbor current positions
% a = forceFunctions(X,x,V,v,10,1,3);
% a_mag = norm(a);

% % construct the current connectivity matrix based on this deformation (only stretched fibers):
% % (NOTE: con_verts_cur is only used by fiber elements)
% con_verts_cur = zeros(num_verts); % current connectivity matrix; symmetrical
% [i_vertsA,i_vertsB] = find(con_verts_ref); % store indices of non-zero elements of ref connectivity matrix
% for i = 1:length(i_vertsA) % iterate through indices of non-zero elements (i.e. check each connection)
%     % calculate edge lengths:
%     edgeLength_ref = sqrt( ...
%         (pos_verts_ref(i_vertsA(i)) - pos_verts_ref(i_vertsB(i)))^2 ...
%         + (pos_verts_ref(i_vertsA(i) + num_verts) - pos_verts_ref(i_vertsB(i)) + num_verts)^2 ...
%         ) / s_f; % scale ref edge length by the inverse of fiber straightness to get the effective ref edge length
%     edgeLength_cur = sqrt( ...
%         (pos_verts_cur(i_vertsA(i)) - pos_verts_cur(i_vertsB(i)))^2 ...
%         + (pos_verts_cur(i_vertsA(i) + num_verts) - pos_verts_cur(i_vertsB(i)) + num_verts)^2 ...
%         );
%     % use Heaviside step function to exclude compressed fibers:
%     con_verts_cur(i_vertsA(i),i_vertsB(i)) = heaviside(edgeLength_cur / edgeLength_ref - 1);
% end



%     % loop once through middle vertices to equilibrate their positions sequentially:
%     F_net_residual = zeros(length(i_verts_group_M),1);
%     for i_vert = 1:length(i_verts_group_M)
%         % vertex of interest:
%         X = [pos_verts_ref(i_verts_group_M(i_vert)), pos_verts_ref(i_verts_group_M(i_vert) + num_verts)]; % vertex's reference position
%         x_guess = [pos_verts_cur(i_verts_group_M(i_vert)), pos_verts_cur(i_verts_group_M(i_vert) + num_verts)]; % starting guess for vertex's current position
%         % neighbor vertices:
%         i_connections = find(con_verts(i_verts_group_M(i_vert),:)); % 1D array containing indices of all the vertex's connections
%         V = [pos_verts_ref(i_connections), pos_verts_ref(i_connections + num_verts)]; % neighbors' reference positions, stored as [vertex, coordinate (x1:x2)]
%         v = [pos_verts_cur(i_connections), pos_verts_cur(i_connections + num_verts)]; % neighbors' current positions, stored as [vertex, coordinate (x1:x2)]
%         % system of equations:
%         F_net_residual(i_vert) = norm(forceFunctions(X,x_guess,V,v,stiffness_f,straightness_f,stiffness_g_t,stiffness_g_c)); % store residual net force (before equilibration)
%         system = @(x)forceFunctions(X,x,V,v,stiffness_f,straightness_f,stiffness_g_t,stiffness_g_c); % functions to solve (see "helper functions" below)
%         x_eq = fsolve(system, x_guess, solver_options); % solve system to find vertex's current equilibrium position
%         % update vertex position:
%         pos_verts_cur(i_verts_group_M(i_vert)) = x_eq(1); % update x1
%         pos_verts_cur(i_verts_group_M(i_vert) + num_verts) = x_eq(2); % update x2
%     end




% F(1) = sum( ...
%     k_f .* ((sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2) - sqrt((V(:,1)-X(1)).^2+(V(:,2)-X(2)).^2) ./ s_f) .* (v(:,1) - x(1)) ./ sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2)) ...
%     .* heaviside(sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2) ./ sqrt((V(:,1)-X(1)).^2+(V(:,2)-X(2)).^2) .* s_f - 1) ...
%     + k_g .* ((sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2) - sqrt((V(:,1)-X(1)).^2+(V(:,2)-X(2)).^2)) .* (v(:,1) - x(1)) ./ sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2)) ...
%     ); % ~ total force on vertex along x1, want to = 0
% F(2) = sum( ...
%     k_f .* ((sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2) - sqrt((V(:,1)-X(1)).^2+(V(:,2)-X(2)).^2) ./ s_f) .* (v(:,2) - x(2)) ./ sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2)) ...
%     .* heaviside(sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2) ./ sqrt((V(:,1)-X(1)).^2+(V(:,2)-X(2)).^2) .* s_f - 1) ...
%     + k_g .* ((sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2) - sqrt((V(:,1)-X(1)).^2+(V(:,2)-X(2)).^2)) .* (v(:,2) - x(2)) ./ sqrt((v(:,1)-x(1)).^2+(v(:,2)-x(2)).^2)) ...
%     ); % ~ total force on vertex along x2, want to = 0


