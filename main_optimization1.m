% Optimized stimulation using multi-electrode for single focal stimulation
clear
rng(200)
num_E = '16'; % # of electrodes
model = 'cylinder'; % 3d model
% 1. cylinder
% 2. sphere : electordes are placed on upper hemisphere
file = ['model_',model,'_',num_E,'.mat'];
load(file) % load the saved the electric field distribution map associated with each model

% Define global constraints and maximum current/power values
I_max = 2; 
P_max = 1;

% Initialize weight arrays for targets and non-target regions
w = zeros(length(data_x),1); 
w_non = zeros(length(data_x),1);

% Define target coordinates and their respective radii
target = [0;0;0]; 
switch model
    case 'cylinder'
        r = 3.5;     
    case 'sphere'
        r = 3;
end

% Calculate spatial weights based on distance to target locations
if length(data_x) <= 10201 +100
    tmp = xyz;
    clear xyz
    xyz.pos = tmp;
    model_bound = 25;
else
    model_bound = (25*0.842);
end
for i = 1:length(w)
    d = xyz.pos(:,i)-target;
    D = norm(d);
    if D < r
        w(i) = (r-D+1)^2;
    elseif norm(xyz.pos(:,i))<= model_bound
        w_non(i) = 1;
    end
end


% Transfer matrix representing the electric field generated per unit current
A = [data_x; data_y; data_z];
A(isnan(A)) = 0;

% lower and upper bounds
lb =  zeros(4*M,1);
I_ub = 1;
ub = ones(4*M,1) * I_ub;

% Initialize random starting points for the optimization solver
a1=rand(M,1); a1= a1*I_ub/sum(a1)/10;
a2=rand(M,1); a2= a2*I_ub/sum(a2)/10;
a3=rand(M,1); a3= a3*I_ub/sum(a3)/10;
a4=rand(M,1); a4= a4*I_ub/sum(a4)/10;
x0 = [a1;a2;a3;a4];

% Linear constraints
a = [ones(1,M) zeros(1,M) ones(1,M) zeros(1,M)];
b = [I_max];
Aeq = [ones(1,M) -ones(1,M) zeros(1,2*M);zeros(1,2*M) ones(1,M) -ones(1,M)];
beq = [0;0];

% optimization option
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1000);
options.Display = 'iter-detailed';
options.MaxFunctionEvaluations = 60000;
nonlcon =@(x)Pmax_constraint(x,A,w_non,P_max,M);

% optimziation
x = fmincon(@(x1)fun(x1,A,w,M),x0,a,b,Aeq,beq,lb,ub,nonlcon,options);
% i_stim : stimulation current [2M x 1]
i_stim = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];
total_I = sum(abs(i_stim))/2;
i_stim = i_stim*2/total_I;
%% Results Visualization
if model == "cylinder"
    E_plot_2d(i_stim,target,r,file);
else
    E_plot_3d(i_stim,target,r,file);
end


%%
% objective function
function out = fun(x,A,w,M)

i_stim = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];
% convert current to electric fields
Ex1 = A(1:length(A)/3,:)*i_stim(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_stim(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_stim(1:M);
E1 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_stim(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_stim(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_stim(M+1:2*M);
E2 = [Ex2 Ey2 Ez2];
% calculate maximum envelope amplitude for each point
E_max = E_MAX(E1,E2);
E_max(isnan(E_max)) = 0;
% objective function to be minimized
out = -w'*((E_max-0).^2);
end

% nonlinear constraints
function [c,ceq] = Pmax_constraint(x,A,w_non,P_max,M)

i_stim = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];
Ex1 = A(1:length(A)/3,:)*i_stim(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_stim(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_stim(1:M);

Ex2 = A(1:length(A)/3,:)*i_stim(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_stim(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_stim(M+1:2*M);

E1 = [Ex1,Ey1,Ez1];
E2 = [Ex2,Ey2,Ez2];
E_max= E_MAX(E1,E2);
E_max(isnan(E_max)) = 0;
% nonlinear constraints
% 1. Power of envelope in non target area < Pmax
c = w_non'*E_max.^2-P_max;
% 2.Currents corresponding to frequencies f and f+df must be separated 
ceq = i_stim(1:M).*i_stim(M+1:2*M);
end

