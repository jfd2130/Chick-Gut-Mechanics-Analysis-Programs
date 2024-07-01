% four required outputs from FEBio simulation (as inidividual files):
% displacement, stress, volume ratio, deformation gradient

% MANUAL INPUTS:
elements_x = 50;
elements_y = 15;
initialLength_x = 2.64;
initialLength_y = 0.82;
elements_endsSkip = 12; % to ignore elements on ends to avoid edge effects of BCs

% define elements we want to keep:
elements_kept = elements_x - 2 * elements_endsSkip;
elementIndices = []; % of the elements we will keep
for i = 1:elements_y
    elementIndices = [elementIndices, 1 + elements_kept*(i-1) + elements_endsSkip*(2*i-1) : elements_kept*i + elements_endsSkip*(2*i-1) ];
end

% import data and extract only those parts we care about:
% displacement:
disp_raw = readmatrix('febioOutput_displacement'); % import node displacement data (x,y,z)
disp_x = max(disp_raw (:, 2:3:end));
disp_y = max(disp_raw (:, 3:3:end));
% stress (Cauchy):
stress_raw  = readmatrix('febioOutput_stress'); % import element stress data (x,y,z,xy,yz,xz)
stress_raw = stress_raw (elementIndices, 2:end); % cut off end elements and element index column
stress = zeros(size(stress_raw, 1), size(stress_raw, 2)/6, 3, 3); % "element, timepoint, stress matrix dim 1, stress matrix dim 2"
for e = 1:size(stress, 1) % element, e
    for t = 1:size(stress, 2) % time, t
        stress(e,t,:,:) =...
            [stress_raw(e,1+6*(t-1)), stress_raw(e,4+6*(t-1)), stress_raw(e,6+6*(t-1));...
            stress_raw(e,4+6*(t-1)), stress_raw(e,2+6*(t-1)), stress_raw(e,5+6*(t-1));...
            stress_raw(e,6+6*(t-1)), stress_raw(e,5+6*(t-1)), stress_raw(e,3+6*(t-1))];
    end
end
% volume ratio (Jacobian, J):
J_raw  = readmatrix('febioOutput_J'); % import element volume ratio (Jacobian) data (J)
J = J_raw(elementIndices, 2:end); % cut off end elements and element index column
% deformation gradient (F):
F_raw  = readmatrix('febioOutput_F'); % import element deformation gradient (F) data (xx,xy,xz,yx,yy,yz,zx,zy,zz);
F_raw = F_raw(elementIndices, 2:end); % cut off end elements and element index column
F = zeros(size(F_raw, 1), size(F_raw, 2)/9, 3, 3); % "element, timepoint, F matrix dim 1, F matrix dim 2"
for e = 1:size(F, 1) % element, e
    for t = 1:size(F, 2) % time, t
        F(e,t,:,:) =...
            [F_raw(e,1+9*(t-1)), F_raw(e,2+9*(t-1)), F_raw(e,3+9*(t-1));...
            F_raw(e,4+9*(t-1)), F_raw(e,5+9*(t-1)), F_raw(e,6+9*(t-1));...
            F_raw(e,7+9*(t-1)), F_raw(e,8+9*(t-1)), F_raw(e,9+9*(t-1))];
    end
end

% calculate 1st Piola-Kirchhoff stress (P):
PK = zeros(size(stress, 1), size(stress, 2), 3, 3);
for e = 1:size(PK, 1) % element, e
    for t = 1:size(PK, 2) % time, t
        PK(e,t,:,:) = squeeze(J(e,t)) * (squeeze(stress(e,t,:,:)) / squeeze(F(e,t,:,:)).' );
    end
end
PK_xx = mean(PK(:,:,1,1));
stress_xx = mean(stress(:,:,1,1));

% take off first timepoint (initial condition step):
initialLength_y = initialLength_y - disp_y(2)*2;
disp_x = disp_x(2:end);
disp_y = disp_y(2:end) - disp_y(2);
PK_xx  = PK_xx (2:end);
stress_xx = stress_xx (2:end);

% compute stretches:
stretch_x = (initialLength_x + disp_x) ./ initialLength_x;
stretch_y = (initialLength_y - disp_y.*2) ./ initialLength_y;

% plot:
figure
plot(stretch_x, PK_xx ,'k','LineWidth',2); hold on
%plot(stretch_x, stress_xx ,'k','LineWidth',1);
plot([1,1.70,2.30,2.43,2.55,2.63,2.70,2.73,2.75],...
    [0,1.28,7.31,12.4,17.4,22.4,27.4,32.5,37.5].*(10^3),...
    ':k','LineWidth',3); % dmso
plot([1,1.34,1.63,1.70,1.76,1.81,1.83,1.84,1.88],...
    [0,1.29,7.39,12.5,17.5,22.5,27.4,32.5,37.5].*(10^3),...
    ':c','LineWidth',3); % cd
plot([1,1.92,2.68,2.88,2.98,3.08,3.18,3.28,3.34],...
    [0,1.33,7.32,12.42,17.46,22.40,27.29,32.40,37.17].*(10^3),...
    ':r','LineWidth',3); % ca
xlabel('axial stretch','FontSize',14);
ylabel('axial stress (Pa)','FontSize',14);
legend({'model', 'dmso', 'cd'});
xlim([1, 4.5]);
ylim([0, 42e3]);

figure
plot(stretch_x, stretch_y,'k','LineWidth',2); hold on
plot([2.83,2.76,2.73,2.39,2.17,1.85,1.89,1.19],...
    [0.333,0.430,0.532,0.624,0.722,0.810,0.904,1.00],...
    ':k','LineWidth',3); % dmso
plot([2.00,1.86,1.80,1.73,1.60,1.43,1.36,1.05],...
    [0.234,0.321,0.433,0.532,0.640,0.745,0.835,0.985],...
    ':c','LineWidth',3); % cd
plot([3.43,3.21,3.05,2.75,2.63,2.45,2.04,1.62],...
    [0.45,0.54,0.62,0.71,0.78,0.87,0.94,1.01],...
    ':r','LineWidth',3); % ca
xlabel('axial stretch','FontSize',14);
ylabel('transverse stretch','FontSize',14);
legend({'model', 'dmso', 'cd'});
xlim([1, 4.5]);
ylim([0, 1.1]);









