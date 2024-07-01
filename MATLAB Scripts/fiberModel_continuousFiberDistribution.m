%% invariable circular distribution:

% inputs:
stretch_axial = 1:0.1:3.2; % axial stretch
stretch_trans = 1 ./ (1 + 0.03259.*(exp((stretch_axial-1).*2.109)-1)); % (DMSO) transverse stretch, calc from fit logistic equation to axial-trans stretch data
%stretch_trans = 1 ./ (1 + 0.06306.*(exp((stretch_axial-1).*3.924)-1)); % (CD, R^2 = )
%stretch_trans = 1 ./ (1 + 0.01197.*(exp((stretch_axial-1).*1.930)-1)); % (CA, R^2 = 0.99)
prestretch = 0.83; % of single fibers (DMSO = 0.83, CD = 0.92)
E_g = 1420; % (ground) modulus in Pa
v_g = 0.19; % (ground) substance Poisson's ratio
d = 100; % number of bins for estimating integral over fibers

% fit values:
E_f = 1200; % (fiber) representing measure of the fiber modulus  <-- FIT
a = 0.1; % (fiber) coefficient of exponential argument  <-- FIT
b = 2; % (fiber) power of exponential argument  <-- FIT

% prestretch:
F_pre = [prestretch 0; 
        0 prestretch];

% predictions and error:
stress = zeros(1,length(stretch_axial));
singleFiberStress = zeros(1,length(stretch_axial));

% predict stress for each axial stretch value:
for i = 1:length(stretch_axial)

    % right Cauchy-Greene tensor (no prestretch):
    F = [stretch_axial(i) 0;
                0 stretch_trans(i)];
    C = [F(1,1)^2 + F(2,1)^2, F(1,2)*F(1,1) + F(2,2)*F(2,1);
        F(1,2)*F(1,1) + F(2,2)*F(2,1), F(2,2)^2 + F(1,2)^2];
    
    % ground substance:
    J = F(1,1)*F(2,2) - F(1,2)*F(2,1); % "volume" ratio
    I_1 = C(1,1) + C(2,2); % first invariant of C
    mu = E_g/(2*(1+v_g)); % a constant of the Neo-Hookean material
    lam = v_g*E_g/((1+v_g)*(1-2*v_g)); % a constant of the Neo-Hookean material
    
    dW_g = mu*F(1,1) + F(2,2)*(lam*log(J)-mu)/J; % dW_g/dF(1,1) ... stress contribution from ground substance

    % right Cauchy-Greene tensor (with prestretch):
    F_applied = [stretch_axial(i) 0;
                0 stretch_trans(i)];
    F = F_applied * F_pre;
    C = [F(1,1)^2 + F(2,1)^2, F(1,2)*F(1,1) + F(2,2)*F(2,1);
        F(1,2)*F(1,1) + F(2,2)*F(2,1), F(2,2)^2 + F(1,2)^2];
    
    % continuous fiber distibution:
    R = 1/pi; % fiber density distribution function (circular distribution)
    I_n = @(theta) cos(theta)^2*C(1,1) + sin(theta)^2*C(2,2) + 2*cos(theta)*sin(theta)*C(1,2);
    dI_n = @(theta) 2 * (cos(theta)^2*F(1,1) + cos(theta)*sin(theta)*F(1,2));
    
    dW_n = @(theta) 2*E_f * exp(a*(I_n(theta)-1)^b) * (I_n(theta)-1)^(b-1) * dI_n(theta); % dW_n/dF(1,1) ... stress contribution from single fiber
    
    dW_f = 0; % dW_f/dF(1,1) ... total stress contribution from fiber distribution; initialized as 0
    for theta = -pi/2:pi/d:pi/2 % ~ approximation of definite integral over [-pi/2,pi/2]
        % calculate at first angle:
        dW_f_1 = 0;
        if I_n(theta) > 1 % ~ Heaviside step function
            dW_f_1 = R * dW_n(theta);
        end
        % calculate at second angle (+pi/d):
        dW_f_2 = 0;
        if theta + pi/d > pi/2 % make sure next step won't take us out of range
            if I_n(pi/2) > 1 % ~ Heaviside step function
                dW_f_1 = R * dW_n(pi/2); % calculate at max range value
            end
        elseif I_n(theta + pi/d) > 1 % ~ Heaviside step function
            dW_f_1 = R * dW_n(theta + pi/d); % calculate at next step value
        end
        % estimate area under curve as rectangle of mean height and pi/d width:
        dW_f = dW_f + mean([dW_f_1 dW_f_2])*(pi/d);
    end
    
    % first Piola-Kirchhoff stress, only P(1,1):
    P_11 = dW_g + dW_f;

    stress(i) = P_11;

end

% visualize:
subplot(1,2,1)
plot(stretch_axial, stretch_trans); hold on
title("input");
xlabel("axial stretch");
ylabel("transverse stretch");
xlim([1,4])
ylim([0,1])
text(2.5, 0.9, {"with λ_{pre} = " string(prestretch)})
subplot(1,2,2)
plot(stretch_axial, stress); hold on
title("output");
xlabel("axial stretch");
ylabel("stress (kPa)");
xlim([1,4])
ylim([1,43e3])

% show data for comparison:
subplot(1,2,2)
plot([1,1.70,2.30,2.43,2.55,2.63,2.70,2.73,2.75],...
    [0,1.28,7.31,12.4,17.4,22.4,27.4,32.5,37.5].*(10^3),...
    ':k','LineWidth',2); hold on % dmso
subplot(1,2,2)
plot([1,1.34,1.63,1.70,1.76,1.81,1.83,1.84,1.88],...
    [0,1.29,7.39,12.5,17.5,22.5,27.4,32.5,37.5].*(10^3),...
    ':', 'Color', [0 0.8 1],'LineWidth',2); hold on % cd


%% variable eliptical distribution:

% inputs:
stretch_axial = 1:0.1:2; % axial stretch
stretch_trans = 1 ./ (1 + 0.03259.*(exp((stretch_axial-1).*2.109)-1)); % (DMSO) transverse stretch, calc from fit logistic equation to axial-trans stretch data
%stretch_trans = 1 ./ (1 + 0.06306.*(exp((stretch_axial-1).*3.924)-1)); % (CD, R^2 = )
%stretch_trans = 1 ./ (1 + 0.01197.*(exp((stretch_axial-1).*1.930)-1)); % (CA, R^2 = 0.99)
prestretch = 0.83; % of single fibers (DMSO = 0.83, CD = 0.92)
E_g = 1420; % (ground) modulus in Pa
v_g = 0.19; % (ground) substance Poisson's ratio
d = 100; % number of bins for estimating integral over fibers

% fit values:
E_f = 1400; % (fiber) representing measure of the fiber modulus  <-- FIT (DMSO ~ 1400)
a = 0.07; % (fiber) coefficient of exponential argument  <-- FIT
b = 2.1; % (fiber) power of exponential argument  <-- FIT

% prestretch:
F_pre = [prestretch 0; 
        0 prestretch];

% predictions and error:
stress = zeros(1,length(stretch_axial));

% predict stress for each axial stretch value:
for i = 1:length(stretch_axial)

    % right Cauchy-Greene tensor (no prestretch):
    F = [stretch_axial(i) 0;
        0 stretch_trans(i)];
    C = [F(1,1)^2 + F(2,1)^2, F(1,2)*F(1,1) + F(2,2)*F(2,1);
        F(1,2)*F(1,1) + F(2,2)*F(2,1), F(2,2)^2 + F(1,2)^2];
    
    % ground substance:
    J = F(1,1)*F(2,2) - F(1,2)*F(2,1); % "volume" ratio
    I_1 = C(1,1) + C(2,2); % first invariant of C
    mu = E_g/(2*(1+v_g)); % a constant of the Neo-Hookean material
    lam = v_g*E_g/((1+v_g)*(1-2*v_g)); % a constant of the Neo-Hookean material
    
    dW_g = mu*F(1,1) + F(2,2)*(lam*log(J)-mu)/J; % dW_g/dF(1,1) ... stress contribution from ground substance

    % right Cauchy-Greene tensor (with prestretch):
    F_applied = [stretch_axial(i) 0;
                0 stretch_trans(i)];
    F = F_applied * F_pre;
    C = [F(1,1)^2 + F(2,1)^2, F(1,2)*F(1,1) + F(2,2)*F(2,1);
        F(1,2)*F(1,1) + F(2,2)*F(2,1), F(2,2)^2 + F(1,2)^2];
    
    % continuous fiber distibution:
    pow = 4.5; % (DMSO ~ 3.4)
    a_R = (F(1,1)/(F(2,2)))^pow; % if b_R = 1, a parameter of R
    da_R = pow * (F(1,1)/(F(2,2)))^(pow-1) * 1/F(2,2);
    e = sqrt(1-(1/a_R)^2); % a parameter of R
    de = da_R / (a_R^3 * e);
    R = @(theta) 1 / (2*ellipke(e) * sqrt((cos(theta)/a_R)^2 + sin(theta)^2)); % fiber density distribution function (eliptical distribution)
    dR = @(theta) (da_R * cos(theta)^2) / (2*ellipke(e) * (a_R*sqrt((cos(theta)/a_R)^2 + sin(theta)^2))^3) - ...
         (ellipticE(e)/(e*(1-e^2)) - ellipke(e)/e) * de / (2*ellipke(e)^2 * sqrt((cos(theta)/a_R)^2 + sin(theta)^2));
    I_n = @(theta) cos(theta)^2*C(1,1) + sin(theta)^2*C(2,2) + 2*cos(theta)*sin(theta)*C(1,2);
    dI_n = @(theta) 2 * (cos(theta)^2*F(1,1) + cos(theta)*sin(theta)*F(1,2));
    W_n = @(theta) E_f/4 * (I_n(theta)-1)^2; % nonlinear fiber:  @(theta) E_f/(a*b) * (exp(a*(I_n(theta)-1)^b)-1);
    dW_n = @(theta) E_f/2 * (I_n(theta)-1) * dI_n(theta); % nonlinear fiber:  @(theta) 2*E_f * exp(a*(I_n(theta)-1)^b) * (I_n(theta)-1)^(b-1) * dI_n(theta); % dW_n/dF(1,1) ... stress contribution from single fiber

    dW_f = 0; % dW_f/dF(1,1) ... total stress contribution from fiber distribution; initialized as 0
    for theta = -pi/2:pi/d:pi/2 % ~ approximation of definite integral over [-pi/2,pi/2]
        % calculate at first angle:
        dW_f_1 = 0;
        if I_n(theta) > 1 % ~ Heaviside step function
            dW_f_1 = dR(theta)*W_n(theta) + R(theta)*dW_n(theta);
        end
        % calculate at second angle (+pi/d):
        dW_f_2 = 0;
        if theta + pi/d > pi/2 % make sure next step won't take us out of range
            if I_n(pi/2) > 1 % ~ Heaviside step function
                dW_f_2 = dR(pi/2)*W_n(pi/2) + R(pi/2)*dW_n(pi/2); % calculate at max range value
            end
        elseif I_n(theta + pi/d) > 1 % ~ Heaviside step function
            dW_f_2 = dR(theta+pi/d)*W_n(theta+pi/d) + R(theta+pi/d)*dW_n(theta+pi/d); % calculate at next step value
        end
        % estimate area under curve as rectangle of mean height and pi/d width:
        dW_f = dW_f + mean([dW_f_1 dW_f_2])*(pi/d);
    end
    
    % first Piola-Kirchhoff stress, but only P(1,1):
    P_11 = dW_g + dW_f;

    stress(i) = P_11;

end

% visualize:
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
subplot(1,2,1)
plot(stretch_axial, stretch_trans); hold on
title("input");
xlabel("axial stretch");
ylabel("transverse stretch");
xlim([1,4])
ylim([0,1])
text(2.5, 0.9, {"with λ_{pre} = " string(prestretch)})
subplot(1,2,2)
plot(stretch_axial, stress); hold on
%plot(((stretch_axial-1).*0.5)+1, stress); % i.e. for model to perfectly predict, we'd need to effectively half all the Δλ... so, working from measured fiber straightnesses, we'd need (0.83/0.92)^7 to produce a result so extreme!
%plot(stretch_axial, stress.*8); % but just scaling stress does not produce the correct curve shape
title("output");
xlabel("axial stretch");
ylabel("stress (kPa)");
xlim([1,4])
ylim([1,43e3])

% show data for comparison:
subplot(1,2,2)
plot([1,1.70,2.30,2.43,2.55,2.63,2.70,2.73,2.75],...
    [0,1.28,7.31,12.4,17.4,22.4,27.4,32.5,37.5].*(10^3),...
    ':k','LineWidth',2); hold on % dmso
subplot(1,2,2)
plot([1,1.34,1.63,1.70,1.76,1.81,1.83,1.84,1.88],...
    [0,1.29,7.39,12.5,17.5,22.5,27.4,32.5,37.5].*(10^3),...
    ':', 'Color', [0 0.8 1],'LineWidth',2); hold on % cd



%% variable eliptical distribution - INITIAL ORTHO:

% inputs:
stretch_axial = 1:0.1:4; % axial stretch
stretch_trans = 1 ./ (1 + 0.03259.*(exp((stretch_axial-1).*2.109)-1)); % (DMSO) transverse stretch, calc from fit logistic equation to axial-trans stretch data
%stretch_trans = 1 ./ (1 + 0.06306.*(exp((stretch_axial-1).*3.924)-1)); % (CD, R^2 = )
%stretch_trans = 1 ./ (1 + 0.01197.*(exp((stretch_axial-1).*1.930)-1)); % (CA, R^2 = 0.99)
prestretch = 0.83; % of single fibers (DMSO = 0.83, CD = 0.92)
E_g = 1420; % (ground) modulus in Pa
v_g = 0.19; % (ground) substance Poisson's ratio
d = 100; % number of bins for estimating integral over fibers

% fit values:
E_f = 270; % (fiber) representing measure of the fiber modulus  <-- FIT
a = 0.1; % (fiber) coefficient of exponential argument  <-- FIT
b = 2.0; % (fiber) power of exponential argument  <-- FIT

% prestretch:
F_pre = [prestretch 0; 
        0 prestretch];

% predictions and error:
stress = zeros(1,length(stretch_axial));

% predict stress for each axial stretch value:
for i = 1:length(stretch_axial)
    % deformation applied for this timestep:
    F_applied = [stretch_axial(i) 0;
                0 stretch_trans(i)];

    % right Cauchy-Greene tensor (no prestretch):
    F = F_applied;
    C = [F(1,1)^2 + F(2,1)^2, F(1,2)*F(1,1) + F(2,2)*F(2,1);
        F(1,2)*F(1,1) + F(2,2)*F(2,1), F(2,2)^2 + F(1,2)^2];
    
    % ground substance:
    J = F(1,1)*F(2,2) - F(1,2)*F(2,1); % "volume" ratio
    I_1 = C(1,1) + C(2,2); % first invariant of C
    mu = E_g/(2*(1+v_g)); % a constant of the Neo-Hookean material
    lam = v_g*E_g/((1+v_g)*(1-2*v_g)); % a constant of the Neo-Hookean material
    
    dW_g = mu*F(1,1) + F(2,2)*(lam*log(J)-mu)/J; % dW_g/dF(1,1) ... stress contribution from ground substance

    % right Cauchy-Greene tensor (with prestretch):
    F = F_applied * F_pre;
    C = [F(1,1)^2 + F(2,1)^2, F(1,2)*F(1,1) + F(2,2)*F(2,1);
        F(1,2)*F(1,1) + F(2,2)*F(2,1), F(2,2)^2 + F(1,2)^2];
    
    % continuous fiber distibution:
    a_R = F(1,1)/(F(2,2)^3); % if b_R = 1, a parameter of R
    da_R = 1/(F(2,2)^3);
    e = sqrt(1-(1/a_R)^2); % a parameter of R ("k" in text)
    de = da_R / (a_R^3 * e);
    R = @(theta) 1 / (2*ellipke(e) * sqrt((cos(theta)/a_R)^2 + sin(theta)^2)); % fiber density distribution function (eliptical distribution)
    dR = @(theta) (da_R * cos(theta)^2) / (2*ellipke(e) * (a_R*sqrt((cos(theta)/a_R)^2 + sin(theta)^2))^3) - ...
         (ellipticE(e)/(e*(1-e^2)) - ellipke(e)/e) * de / (2*ellipke(e)^2 * sqrt((cos(theta)/a_R)^2 + sin(theta)^2));
    I_n = @(theta) cos(theta)^2*C(1,1) + sin(theta)^2*C(2,2) + 2*cos(theta)*sin(theta)*C(1,2);
    dI_n = @(theta) 2 * (cos(theta)^2*F(1,1) + cos(theta)*sin(theta)*F(1,2));
    W_n = @(theta) E_f/(a*b) * (exp(a*(I_n(theta)-1)^b)-1); % linear: @(theta) E_f/4 * (I_n(theta)-1)^2;
    dW_n = @(theta) 2*E_f * exp(a*(I_n(theta)-1)^b) * (I_n(theta)-1)^(b-1) * dI_n(theta); % linear: @(theta) E_f/2 * (I_n(theta)-1) * dI_n(theta); % dW_n/dF(1,1) ... stress contribution from single fiber

    dW_f = 0; % dW_f/dF(1,1) ... total stress contribution from fiber distribution; initialized as 0
    for theta = -pi/2:pi/d:pi/2 % ~ approximation of definite integral over [-pi/2,pi/2]
        % calculate at first angle:
        dW_f_1 = 0;
        if I_n(theta) > 1 % ~ Heaviside step function
            dW_f_1 = dR(theta)*W_n(theta) + R(theta)*dW_n(theta);
        end
        % calculate at second angle (+pi/d):
        dW_f_2 = 0;
        if theta + pi/d > pi/2 % make sure next step won't take us out of range
            if I_n(pi/2) > 1 % ~ Heaviside step function
                dW_f_2 = dR(pi/2)*W_n(pi/2) + R(pi/2)*dW_n(pi/2); % calculate at max range value
            end
        elseif I_n(theta + pi/d) > 1 % ~ Heaviside step function
            dW_f_2 = dR(theta+pi/d)*W_n(theta+pi/d) + R(theta+pi/d)*dW_n(theta+pi/d); % calculate at next step value
        end
        % estimate area under curve as rectangle of mean height and pi/d width:
        dW_f = dW_f + mean([dW_f_1 dW_f_2])*(pi/d);
    end
    
    % first Piola-Kirchhoff stress, but only P(1,1):
    P_11 = dW_g + dW_f;

    stress(i) = P_11;

end

% visualize:
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
subplot(1,2,1)
plot(stretch_axial, stretch_trans); hold on
title("input");
xlabel("axial stretch");
ylabel("transverse stretch");
xlim([1,4])
ylim([0,1])
text(2.5, 0.9, {"with λ_{pre} = " string(prestretch)})
subplot(1,2,2)
plot(stretch_axial, stress); hold on
%plot(((stretch_axial-1).*0.5)+1, stress); % i.e. for model to perfectly predict, we'd need to effectively half all the Δλ... so, working from measured fiber straightnesses, we'd need (0.83/0.92)^7 to produce a result so extreme!
%plot(stretch_axial, stress.*8); % but just scaling stress does not produce the correct curve shape
title("output");
xlabel("axial stretch");
ylabel("stress (kPa)");
xlim([1,4])
ylim([1,43e3])

% show data for comparison:
subplot(1,2,2)
plot([1,1.70,2.30,2.43,2.55,2.63,2.70,2.73,2.75],...
    [0,1.28,7.31,12.4,17.4,22.4,27.4,32.5,37.5].*(10^3),...
    ':k','LineWidth',2); hold on % dmso
subplot(1,2,2)
plot([1,1.34,1.63,1.70,1.76,1.81,1.83,1.84,1.88],...
    [0,1.29,7.39,12.5,17.5,22.5,27.4,32.5,37.5].*(10^3),...
    ':', 'Color', [0 0.8 1],'LineWidth',2); hold on % cd


%% invariable von Mises distribution:

% inputs:
stretch_axial = 1:0.1:4.5; % axial stretch
stretch_trans = 1 ./ (1 + 0.03259.*(exp((stretch_axial-1).*2.109)-1)); % (DMSO) transverse stretch, calc from fit logistic equation to axial-trans stretch data
%stretch_trans = 1 ./ (1 + 0.06306.*(exp((stretch_axial-1).*3.924)-1)); % (CD, R^2 = )
%stretch_trans = 1 ./ (1 + 0.01197.*(exp((stretch_axial-1).*1.930)-1)); % (CA, R^2 = 0.99)
prestretch = 0.83; % of single fibers (DMSO = 0.83, CD = 0.92)
E_g = 1420; % (ground) modulus in Pa
v_g = 0.19; % (ground) substance Poisson's ratio
d = 100; % number of bins for estimating integral over fibers

% fit values:
E_f = 500000000000000; % (fiber) representing measure of the fiber modulus  <-- FIT (DMSO ~ 1400)
a = 0.05; % (fiber) coefficient of exponential argument  <-- FIT
b = 2; % (fiber) power of exponential argument  <-- FIT

% prestretch:
F_pre = [prestretch 0; 
        0 prestretch];

% predictions and error:
stress = zeros(1,length(stretch_axial));

% predict stress for each axial stretch value:
for i = 1:length(stretch_axial)

    % right Cauchy-Greene tensor (no prestretch):
    F = [stretch_axial(i) 0;
        0 stretch_trans(i)];
    C = [F(1,1)^2 + F(2,1)^2, F(1,2)*F(1,1) + F(2,2)*F(2,1);
        F(1,2)*F(1,1) + F(2,2)*F(2,1), F(2,2)^2 + F(1,2)^2];
    
    % ground substance:
    J = F(1,1)*F(2,2) - F(1,2)*F(2,1); % "volume" ratio
    I_1 = C(1,1) + C(2,2); % first invariant of C
    mu = E_g/(2*(1+v_g)); % a constant of the Neo-Hookean material
    lam = v_g*E_g/((1+v_g)*(1-2*v_g)); % a constant of the Neo-Hookean material
    
    dW_g = mu*F(1,1) + F(2,2)*(lam*log(J)-mu)/J; % dW_g/dF(1,1) ... stress contribution from ground substance

    % right Cauchy-Greene tensor (with prestretch):
    F_applied = [stretch_axial(i) 0;
                0 stretch_trans(i)];
    F_f = F_applied * F_pre;
    C = [F_f(1,1)^2 + F_f(2,1)^2, F_f(1,2)*F_f(1,1) + F_f(2,2)*F_f(2,1);
        F_f(1,2)*F_f(1,1) + F_f(2,2)*F_f(2,1), F_f(2,2)^2 + F_f(1,2)^2];
    
    % continuous fiber distibution:
    b_R = 3; % concentration parameter for 2D von Mises distribution, b_R>0 (b_R=0 is circular)
    bessel_0 = besseli(0, b_R); % modified Bessel function of the first kind of order 0 for b_R
    R = @(theta) exp(b_R * (2*sin(theta)^2 - 1)) / (pi * bessel_0); % fiber density distribution function (2D von Mises distribution, along n_2)
    dR = @(theta) 0;
    I_n = @(theta) cos(theta)^2*C(1,1) + sin(theta)^2*C(2,2) + 2*cos(theta)*sin(theta)*C(1,2);
    dI_n = @(theta) 2 * (cos(theta)^2*F_f(1,1) + cos(theta)*sin(theta)*F_f(1,2));
    W_n = @(theta) E_f/4 * (I_n(theta)-1)^2; % nonlinear fiber:  @(theta) E_f/(a*b) * (exp(a*(I_n(theta)-1)^b)-1);
    dW_n = @(theta) E_f/2 * (I_n(theta)-1) * dI_n(theta); % nonlinear fiber:  @(theta) 2*E_f * exp(a*(I_n(theta)-1)^b) * (I_n(theta)-1)^(b-1) * dI_n(theta); % dW_n/dF(1,1) ... stress contribution from single fiber

    dW_f = 0; % dW_f/dF(1,1) ... total stress contribution from fiber distribution; initialized as 0
    for theta = -pi/2:pi/d:pi/2 % ~ approximation of definite integral over [-pi/2,pi/2]
        % calculate at first angle:
        dW_f_1 = 0;
        if I_n(theta) > 1 % ~ Heaviside step function
            dW_f_1 = dR(theta)*W_n(theta) + R(theta)*dW_n(theta);
        end
        % calculate at second angle (+pi/d):
        dW_f_2 = 0;
        if theta + pi/d > pi/2 % make sure next step won't take us out of range
            if I_n(pi/2) > 1 % ~ Heaviside step function
                dW_f_2 = dR(pi/2)*W_n(pi/2) + R(pi/2)*dW_n(pi/2); % calculate at max range value
            end
        elseif I_n(theta + pi/d) > 1 % ~ Heaviside step function
            dW_f_2 = dR(theta+pi/d)*W_n(theta+pi/d) + R(theta+pi/d)*dW_n(theta+pi/d); % calculate at next step value
        end
        % estimate area under curve as rectangle of mean height and pi/d width:
        dW_f = dW_f + mean([dW_f_1 dW_f_2])*(pi/d);
    end
    
    % first Piola-Kirchhoff stress, but only P(1,1):
    P_11 = dW_g + dW_f;

    stress(i) = P_11;

end

% visualize:
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
subplot(1,2,1)
plot(stretch_axial, stretch_trans); hold on
title("input");
xlabel("axial stretch");
ylabel("transverse stretch");
xlim([1,4])
ylim([0,1])
text(2.5, 0.9, {"with λ_{pre} = " string(prestretch)})
subplot(1,2,2)
plot(stretch_axial, stress); hold on
%plot(((stretch_axial-1).*0.5)+1, stress); % i.e. for model to perfectly predict, we'd need to effectively half all the Δλ... so, working from measured fiber straightnesses, we'd need (0.83/0.92)^7 to produce a result so extreme!
%plot(stretch_axial, stress.*8); % but just scaling stress does not produce the correct curve shape
title("output");
xlabel("axial stretch");
ylabel("stress (kPa)");
xlim([1,4])
ylim([1,43e3])

% show data for comparison:
subplot(1,2,2)
plot([1,1.70,2.30,2.43,2.55,2.63,2.70,2.73,2.75],...
    [0,1.28,7.31,12.4,17.4,22.4,27.4,32.5,37.5].*(10^3),...
    ':k','LineWidth',2); hold on % dmso
subplot(1,2,2)
plot([1,1.34,1.63,1.70,1.76,1.81,1.83,1.84,1.88],...
    [0,1.29,7.39,12.5,17.5,22.5,27.4,32.5,37.5].*(10^3),...
    ':', 'Color', [0 0.8 1],'LineWidth',2); hold on % cd


%% VARYING von Mises distribution:

% inputs:
condition = 1; % 0=DMSO, 1=CD, 2=CA
stretch_axial = 1:0.1:4; % axial stretch
b_base = 0.5; % initial concentration parameter for von Mises distribution, at F(1,1)/F(2,2)=1
defRatioAtCirc = 2.2; % F(1,1)/F(2,2) ratio (>=1) at which von Mises distribution reaches a circle (i.e. concentration parameter = 0), = 1 means starts as circle
E_g = 1420; % (ground) modulus in Pa
v_g = 0.19; % (ground) substance Poisson's ratio
d = 100; % number of bins for estimating integral over fibers

% fit values:
E_f = 3500; % (fiber) representing measure of the fiber modulus  <-- FIT (DMSO ~ 1400)
a = 0.1; % (fiber) coefficient of exponential argument  <-- FIT
b = 2; % (fiber) power of exponential argument  <-- FIT

% initialize condition-specific parameters:
if condition == 0
    stretch_trans = 1 ./ (1 + 0.03259.*(exp((stretch_axial-1).*2.109)-1)); % (DMSO, R^2 = 0.95) transverse stretch, calc from fit equation to axial-trans stretch data
    prestretch = 0.83;
elseif condition == 1
    stretch_trans = 1 ./ (1 + 0.06306.*(exp((stretch_axial-1).*3.924)-1)); % (CD, R^2 = 0.99)
    prestretch = 0.92;
elseif condition == 2
    stretch_trans = 1 ./ (1 + 0.01197.*(exp((stretch_axial-1).*1.930)-1)); % (CA, R^2 = 0.99)
    prestretch = 0.0; % <-- UNKOWN
else
    disp("ERROR: no such condition")
end

% prestretch:
F_pre = [prestretch 0;
        0 prestretch];

% predictions and error:
stress = zeros(1,length(stretch_axial));

% predict stress for each axial stretch value:
for i = 1:length(stretch_axial)
    % ---
    % deformation applied for this timestep:
    F = [stretch_axial(i) 0;
        0 stretch_trans(i)];

    % ---
    % right Cauchy-Greene tensor (no prestretch):
    F_g = F;
    C = [F_g(1,1)^2 + F_g(2,1)^2, F_g(1,2)*F_g(1,1) + F_g(2,2)*F_g(2,1);
        F_g(1,2)*F_g(1,1) + F_g(2,2)*F_g(2,1), F_g(2,2)^2 + F_g(1,2)^2];
    
    % ground substance:
    J = F_g(1,1)*F_g(2,2) - F_g(1,2)*F_g(2,1); % "volume" ratio
    I_1 = C(1,1) + C(2,2); % first invariant of C
    mu = E_g/(2*(1+v_g)); % a constant of the Neo-Hookean material
    lam = v_g*E_g/((1+v_g)*(1-2*v_g)); % a constant of the Neo-Hookean material
    
    dW_g = mu*F_g(1,1) + F_g(2,2)*(lam*log(J)-mu)/J; % dW_g/dF(1,1) ... stress contribution from ground substance

    % ---
    % right Cauchy-Greene tensor (with prestretch):
    F_f = F * F_pre;
    C = [F_f(1,1)^2 + F_f(2,1)^2, F_f(1,2)*F_f(1,1) + F_f(2,2)*F_f(2,1);
        F_f(1,2)*F_f(1,1) + F_f(2,2)*F_f(2,1), F_f(2,2)^2 + F_f(1,2)^2];
    
    % continuous fiber distibution:
    if F(1,1)/F(2,2) <= defRatioAtCirc
        b_R = b_base * (defRatioAtCirc - F(1,1)/F(2,2))/(defRatioAtCirc - 1); % concentration parameter for 2D von Mises distribution, b_R>0 (b_R=0 is circular)
        db_R = -b_base/F(2,2);
        bsl_0 = besseli(0, b_R); % modified Bessel function of the first kind of order 0 for b_R, I_0(b_R)
        bsl_1 = besseli(1, b_R); % modified Bessel function of the first kind of order 1 for b_R, I_1(b_R)
        R = @(theta) exp(b_R * (2*sin(theta)^2-1)) / (pi * bsl_0); % fiber density distribution function (2D von Mises distribution, along n_2)
        dR = @(theta) exp(b_R * (2*sin(theta)^2-1)) / (pi * bsl_0^2) * ((2*sin(theta)^2-1)*bsl_0 - bsl_1) * db_R;
    else
        b_R = b_base * (F(1,1)/F(2,2) - defRatioAtCirc)/(defRatioAtCirc - 1); % concentration parameter for 2D von Mises distribution, b_R>0 (b_R=0 is circular)
        db_R = b_base/F(2,2);
        bsl_0 = besseli(0, b_R); % modified Bessel function of the first kind of order 0 for b_R, I_0(b_R)
        bsl_1 = besseli(1, b_R); % modified Bessel function of the first kind of order 1 for b_R, I_1(b_R)
        R = @(theta) exp(b_R * (2*cos(theta)^2-1)) / (pi * bsl_0); % fiber density distribution function (2D von Mises distribution, along n_2)
        dR = @(theta) exp(b_R * (2*cos(theta)^2-1)) / (pi * bsl_0^2) * ((2*cos(theta)^2-1)*bsl_0 - bsl_1) * db_R;
    end
    I_n = @(theta) cos(theta)^2*C(1,1) + sin(theta)^2*C(2,2) + 2*cos(theta)*sin(theta)*C(1,2);
    dI_n = @(theta) 2 * (cos(theta)^2*F_f(1,1) + cos(theta)*sin(theta)*F_f(1,2));
    W_n = @(theta) E_f/4 * (I_n(theta)-1)^2; % linear fiber
    dW_n = @(theta) E_f/2 * (I_n(theta)-1) * dI_n(theta); % linear fiber; dW_n/dF(1,1) ... stress contribution from single fiber
    %W_n = @(theta) E_f/(a*b) * (exp(a*(I_n(theta)-1)^b)-1); % nonlinear fiber
    %dW_n = @(theta) 2*E_f * exp(a*(I_n(theta)-1)^b) * (I_n(theta)-1)^(b-1) * dI_n(theta); % nonlinear fiber

    dW_f = 0; % dW_f/dF(1,1) ... total stress contribution from fiber distribution; initialized as 0
    for theta = -pi/2:pi/d:pi/2 % ~ approximation of definite integral over [-pi/2,pi/2]
        % calculate at first angle:
        dW_f_1 = 0;
        if I_n(theta) > 1 % ~ Heaviside step function
            dW_f_1 = dR(theta)*W_n(theta) + R(theta)*dW_n(theta);
        end
        % calculate at second angle (+pi/d):
        dW_f_2 = 0;
        if theta + pi/d > pi/2 % make sure next step won't take us out of range
            if I_n(pi/2) > 1 % ~ Heaviside step function
                dW_f_2 = dR(pi/2)*W_n(pi/2) + R(pi/2)*dW_n(pi/2); % calculate at max range value
            end
        elseif I_n(theta + pi/d) > 1 % ~ Heaviside step function
            dW_f_2 = dR(theta+pi/d)*W_n(theta+pi/d) + R(theta+pi/d)*dW_n(theta+pi/d); % calculate at next step value
        end
        % estimate area under curve as rectangle of mean height and pi/d width:
        dW_f = dW_f + mean([dW_f_1 dW_f_2])*(pi/d);
    end
    
    % ---
    % first Piola-Kirchhoff stress, but only P(1,1):
    P_11 = dW_g + dW_f;

    stress(i) = P_11;

end

% visualize:
figure
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
subplot(1,2,1)
plot(stretch_axial, stretch_trans, 'Color', [0 0 0]); hold on
title("input");
xlabel("axial stretch");
ylabel("transverse stretch");
xlim([1,4])
ylim([0,1])
text(2.5, 0.9, {"with λ_{pre} = " string(prestretch)})
subplot(1,2,2)
plot(stretch_axial, stress, 'Color', [0 0 0]); hold on
%plot(((stretch_axial-1).*0.5)+1, stress); % i.e. for model to perfectly predict, we'd need to effectively half all the Δλ... so, working from measured fiber straightnesses, we'd need (0.83/0.92)^7 to produce a result so extreme!
%plot(stretch_axial, stress.*8); % but just scaling stress does not produce the correct curve shape
title("output");
xlabel("axial stretch");
ylabel("stress (kPa)");
xlim([1,4])
ylim([1,43e3])

% show data for comparison:
subplot(1,2,2)
plot([1,1.70,2.30,2.43,2.55,2.63,2.70,2.73,2.75],...
    [0,1.28,7.31,12.4,17.4,22.4,27.4,32.5,37.5].*(10^3),...
    ':k','LineWidth',2); hold on % dmso
subplot(1,2,2)
plot([1,1.34,1.63,1.70,1.76,1.81,1.83,1.84,1.88],...
    [0,1.29,7.39,12.5,17.5,22.5,27.4,32.5,37.5].*(10^3),...
    ':', 'Color', [0 0.8 1],'LineWidth',2); hold on % cd


%% VARYING von Mises distribution (with kappa equation from SHG data):

% inputs:
condition = 1; % 0=DMSO, 1=CD, 2=CA
stretch_axial = 1:0.1:4; % axial stretch
E_g = 1420; % (ground) modulus in Pa, of DMSO toe
v_g = 0.19; % (ground) substance Poisson's ratio, of DMSO toe
d = 100; % number of bins for estimating integral over fibers
% fit values:
E_f = 300; % (fiber) representing measure of the fiber modulus  <-- FIT (DMSO ~ 1400)
a = 0.1; % (fiber) coefficient of exponential argument  <-- FIT
b = 2; % (fiber) power of exponential argument  <-- FIT

% initialize condition-specific parameters:
if condition == 0
    stretch_trans = 1 ./ (1 + 0.03259.*(exp((stretch_axial-1).*2.109)-1)); % (DMSO, R^2 = 0.95) transverse stretch, calc from fit logistic equation to axial-trans stretch data
    prestretch = 0.83;
elseif condition == 1
    stretch_trans = 1 ./ (1 + 0.06306.*(exp((stretch_axial-1).*3.924)-1)); % (CD, R^2 = 0.99)
    prestretch = 0.92;
elseif condition == 2
    stretch_trans = 1 ./ (1 + 0.01197.*(exp((stretch_axial-1).*1.930)-1)); % (CA, R^2 = 0.99)
    prestretch = 0.0; % <-- UNKOWN
else
    disp("ERROR: no such condition")
end

% prestretch deformation matrix:
F_pre = [prestretch 0;
        0 prestretch];

% initialize predictions and error:
stress = zeros(1,length(stretch_axial));

% predict stress for each axial stretch value:
for i = 1:length(stretch_axial)
    % deformation applied for this timestep:
    F = [stretch_axial(i) 0;
        0 stretch_trans(i)];

    % GROUND:
    % right Cauchy-Greene tensor (no prestretch):
    F_g = F;
    C = [F_g(1,1)^2 + F_g(2,1)^2, F_g(1,2)*F_g(1,1) + F_g(2,2)*F_g(2,1);
        F_g(1,2)*F_g(1,1) + F_g(2,2)*F_g(2,1), F_g(2,2)^2 + F_g(1,2)^2];
    
    % ground substance:
    J = F_g(1,1)*F_g(2,2) - F_g(1,2)*F_g(2,1); % "volume" ratio
    I_1 = C(1,1) + C(2,2); % first invariant of C
    mu = E_g/(2*(1+v_g)); % a constant of the Neo-Hookean material
    lam = v_g*E_g/((1+v_g)*(1-2*v_g)); % a constant of the Neo-Hookean material
    
    dW_g = mu*F_g(1,1) + F_g(2,2)*(lam*log(J)-mu)/J; % dW_g/dF(1,1) ... stress contribution from ground substance

    % FIBERS:
    % right Cauchy-Greene tensor (with prestretch):
    F_f = F * F_pre;
    C = [F_f(1,1)^2 + F_f(2,1)^2, F_f(1,2)*F_f(1,1) + F_f(2,2)*F_f(2,1);
        F_f(1,2)*F_f(1,1) + F_f(2,2)*F_f(2,1), F_f(2,2)^2 + F_f(1,2)^2];
    
    % continuous fiber distibution:
    %k = 6.101 / (1 + exp(-F(1,1)/F(2,2))) - 4.82; % logistic function for von Mises concentration parameter, kappa (fit to SHG alignment data for DMSO vs F11/F22)
    %dk = 6.101 * exp(-F(1,1)/F(2,2)) / (F(2,2) * (1+exp(-F(1,1)/F(2,2)))^2);
    k = 7.735 / (1 + exp(-F(1,1))) - 6.036; % logistic function for von Mises concentration parameter, kappa (fit to SHG alignment data for DMSO & CD vs F11)
    dk = 7.735 * exp(-F(1,1)) / (1+exp(-F(1,1)))^2;
    R = @(theta) exp(k*cos(2*theta)) / (pi*besseli(0,k));
    dR = @(theta) exp(k*cos(2*theta)) / (pi*besseli(0,k)) * (cos(2*theta)-besseli(1,k)/besseli(0,k)) * dk;
    I_n = @(theta) cos(theta)^2*C(1,1) + sin(theta)^2*C(2,2) + 2*cos(theta)*sin(theta)*C(1,2);
    dI_n = @(theta) 2 * (cos(theta)^2*F_f(1,1) + cos(theta)*sin(theta)*F_f(1,2));
    %W_n = @(theta) E_f/4 * (I_n(theta)-1)^2; % linear fiber
    %dW_n = @(theta) E_f/2 * (I_n(theta)-1) * dI_n(theta); % linear fiber; dW_n/dF(1,1) ... stress contribution from single fiber
    W_n = @(theta) E_f/(a*b) * (exp(a*(I_n(theta)-1)^b)-1); % nonlinear fiber
    dW_n = @(theta) 2*E_f * exp(a*(I_n(theta)-1)^b) * (I_n(theta)-1)^(b-1) * dI_n(theta); % nonlinear fiber
    
%     figure(1)
    %scatter(F(1,1), dk); hold on
    %figure(2)
    %plot(-pi/2:pi/d:pi/2, R(-pi/2:pi/d:pi/2), 'Color', [i/length(stretch_axial) 0.5 0.5]); hold on
%     W_plot = [];
%     for theta = -pi/2:pi/d:pi/2
%            W_plot(end+1) = dR(theta);
%     end
%     plot3(-pi/2:pi/d:pi/2, zeros(length(W_plot))+i, W_plot,  'Color', [i/length(stretch_axial) 0.5 0.5]); hold on

    dW_f = 0; % dW_f/dF(1,1) ... total stress contribution from fiber distribution; initialized as 0
    for theta = -pi/2:pi/d:pi/2 % ~ approximation of definite integral over [-pi/2,pi/2]
        % calculate at first angle:
        dW_f_1 = 0;
        if I_n(theta) > 1 % ~ Heaviside step function
            dW_f_1 = dR(theta)*W_n(theta) + R(theta)*dW_n(theta);
        end
        % calculate at second angle (+pi/d):
        dW_f_2 = 0;
        if theta + pi/d > pi/2 % make sure next step won't take us out of range
            if I_n(pi/2) > 1 % ~ Heaviside step function
                dW_f_2 = dR(pi/2)*W_n(pi/2) + R(pi/2)*dW_n(pi/2); % calculate at max range value
            end
        elseif I_n(theta + pi/d) > 1 % ~ Heaviside step function
            dW_f_2 = dR(theta+pi/d)*W_n(theta+pi/d) + R(theta+pi/d)*dW_n(theta+pi/d); % calculate at next step value
        end
        % estimate area under curve as rectangle of mean height and pi/d width:
        dW_f = dW_f + mean([dW_f_1 dW_f_2])*(pi/d);
    end
    
    % TOTAL:
    % first Piola-Kirchhoff stress, but only want P(1,1):
    P_11 = dW_g + dW_f;
    stress(i) = P_11;

end

% visualize:
figure
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
subplot(1,2,1)
plot(stretch_axial, stretch_trans, 'Color', [0 0 0]); hold on
title("input");
xlabel("axial stretch");
ylabel("transverse stretch");
xlim([1,4])
ylim([0,1])
text(2.5, 0.9, {"with λ_{pre} = " string(prestretch)})
subplot(1,2,2)
plot(stretch_axial, stress, 'Color', [0 0 0]); hold on
%plot(((stretch_axial-1).*0.5)+1, stress); % i.e. for model to perfectly predict, we'd need to effectively half all the Δλ... so, working from measured fiber straightnesses, we'd need (0.83/0.92)^7 to produce a result so extreme!
%plot(stretch_axial, stress.*8); % but just scaling stress does not produce the correct curve shape
title("output");
xlabel("axial stretch");
ylabel("stress (kPa)");
xlim([1,4])
ylim([1,43e3])

% show data for comparison:
subplot(1,2,2)
plot([1,1.70,2.30,2.43,2.55,2.63,2.70,2.73,2.75],...
    [0,1.28,7.31,12.4,17.4,22.4,27.4,32.5,37.5].*(10^3),...
    ':k','LineWidth',2); hold on % dmso
subplot(1,2,2)
plot([1,1.34,1.63,1.70,1.76,1.81,1.83,1.84,1.88],...
    [0,1.29,7.39,12.5,17.5,22.5,27.4,32.5,37.5].*(10^3),...
    ':', 'Color', [0 0.8 1],'LineWidth',2); hold on % cd








