clear; clc; close all;
addpath('./QuadraticElement');
load('beamsData.mat');

% min m1 + m2 + m3 + m4 + m5
% s.t. deflx < maxDeflx
%      defly < maxDefly
%  u = inv(K)*f

%% General data
d = 2000; % Section width - [mm]
t = 5;  % Section thickness - [mm]

E = 210e3;  % [MPa]
rho = 7850; % [kg/m^3]
A = 4*d*t; % = ((t+d).^2 - (t-d).^2)         % [mm^2]
J = 2/3*d^3*t; % ~ ((t+d).^4 - (t-d).^4)/12; % [mm^4]

% Beam initialization
beams{1} = MyBeams(E,rho,A,J,1000*b0.len{1},b0.ang{1}(2));
beams{2} = MyBeams(E,rho,A,J,1000*b0.len{2},b0.ang{2}(2));
beams{3} = MyBeams(E,rho,A,J,1000*b0.len{3},b0.ang{3}(2));
beams{4} = MyBeams(E,rho,A,J,1000*b0.len{4},b0.ang{4}(2));
beams{5} = MyBeams(E,rho,A,J,1000*b0.len{5},b0.ang{5}(2));

nBeams = numel(beams);

% Constraint initialization
constr{1} = MyConstr(nBeams,"hinge",1,1);     % Node B
constr{2} = MyConstr(nBeams,"hinge",3,1);     % Node A
constr{3} = MyConstr(nBeams,"clamp",1,3,2,1); % Node F' & F''
constr{4} = MyConstr(nBeams,"hinge",3,3,4,1); % Node E' & E''
constr{5} = MyConstr(nBeams,"clamp",4,3,5,1); % Node C' & C'''
constr{6} = MyConstr(nBeams,"hinge",2,3,5,1); % Node C'' & C'''
constr{7} = MyConstr(nBeams,"hinge",1,3);     % Node F'

% Forces initialization
forces{1} = MyForces(nBeams,5,3,0,-30e3*9.8086,0); % On beam 5, 3rd node, put a Fx = 0, Fy = -30e3*9.8086, Mz = 0;

%%
FEM = MyFEM(beams,constr,forces);

%% Node D displacement
[uD,wD,thD] = FEM.getNodeDisplacement(5,2);

defl = sqrt(uD.^2 + wD.^2);
% %% ? vincular reaction on G
% FxG = full(FEM.l(17));
% FyG = full(FEM.l(18));
% 
% FG = beams{6}.R * [FxG; FyG; 0];