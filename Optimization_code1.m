% Optimized stimulation using multi-electrode
clear
num_E = '32'; % # of electrodes
model = 'petridish'; % 3d model
% 1. petridish
% 2. sphere_xy : electordes are placed on the equator of the sphere
% 3. sphere_uphemi : electordes are placed on upper hemisphere
file = ['save_',model,'_',num_E,'.mat'];
load(file) % load the saved the electric field distribution map associated with each model
rng(200)

I_max = 2; % or 1; % maximum current intensity
P_max = 1; % Power limitation in non target area
e = zeros(length(data_x),1); % weight parameter for target area
gamma = zeros(length(data_x),1);% weight parameter for non target area

target = [0;0;0]; % target point
switch model
    case 'petridish'
        r = 3.5; % radius of target area
    case 'sphere_uphemi'
        r = 3;
end
v = [1;0;0];
v = v/norm(v); % a vector of interest, use when want to know the amplitude in a particular direction

%%%%% Set weight for target and non target
if length(data_x) <= 10201 +100
    tmp = xyz;
    clear xyz
    xyz.pos = tmp;
    model_bound = 25;
else
    model_bound = (25*0.842);
end
for i = 1:length(e)
    d = xyz.pos(:,i)-target;
    D = norm(d);
    if D < r
        e(i) = (r-D+1)^2;
    elseif norm(xyz.pos(:,i))<= model_bound
        gamma(i) = 1;
    end
end
%%%%%%

% Electric field distribution map generated per unit current at one electrode
A=[data_x;data_y;data_z];
A(isnan(A)) = 0;

% lower bound and upper bound
lb =  zeros(4*M,1);
I_ub = 1;
ub = ones(4*M,1) * I_ub;

% random initial value
a1=rand(M,1); a1= a1*I_ub/sum(a1)/10;
a2=rand(M,1); a2= a2*I_ub/sum(a2)/10;
a3=rand(M,1); a3= a3*I_ub/sum(a3)/10;
a4=rand(M,1); a4= a4*I_ub/sum(a4)/10;
x0 = [a1;a2;a3;a4];

% Linear constraints
a = [ones(1,M) zeros(1,M) ones(1,M) zeros(1,M)];
b = [I_max]; % a*x < b
Aeq = [ones(1,M) -ones(1,M) zeros(1,2*M);zeros(1,2*M) ones(1,M) -ones(1,M)];
beq = [0;0]; % Aeq*x = beq
%%
% optimization option
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1000);
options.Display = 'iter-detailed';
options.MaxFunctionEvaluations = 60000;
nonlcon =@(x)Pmax_constraint(x,A,gamma,P_max,M);
% optimziation
x = fmincon(@(x1)fun(x1,A,e,gamma,M),x0,a,b,Aeq,beq,lb,ub,nonlcon,options);
% s : optimized current [2M x 1], s(1~M) is current of f, s(M+1~2M) is current of f+df.  
s = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];

%%
total_I = sum(abs(s))/2;
s = s*2/total_I;
% plot envelope amplitude distribution
if length(data_x) <= 10201 +100
    E_plot_2d_flip(s,target,r,v,4,file);
else
    E_plot(s,target,r,v,1,file);
end

% save_file = ['result_opti_code1_3d_000.mat'];
% save(save_file)
%%
% objective function
function out = fun(x,A,e,gamma,M)

s = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];
% convert current to electric fields
Ex1 = A(1:length(A)/3,:)*s(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*s(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*s(1:M);
E1 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*s(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*s(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*s(M+1:2*M);
E2 = [Ex2 Ey2 Ez2];
% calculate maximum envelope amplitude for each point
E_max = E_MAX(E1,E2);
E_max(isnan(E_max)) = 0;
% objective function to be minimized
out = -e'*((E_max-0).^2);
end

% nonlinear constraints
function [c,ceq] = Pmax_constraint(x,A,gamma,P_max,M)

s = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];
Ex1 = A(1:length(A)/3,:)*s(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*s(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*s(1:M);

Ex2 = A(1:length(A)/3,:)*s(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*s(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*s(M+1:2*M);

E1 = [Ex1,Ey1,Ez1];
E2 = [Ex2,Ey2,Ez2];
E_max= E_MAX(E1,E2);
E_max(isnan(E_max)) = 0;
% nonlinear constraints
% 1. Power of envelope in non target area < Pmax
c = gamma'*E_max.^2-P_max;
% 2.Currents corresponding to frequencies f and f+df must be separated 
ceq = s(1:M).*s(M+1:2*M);
end

