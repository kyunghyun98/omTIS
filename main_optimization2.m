% Optimized stimulation using multi-electrode for multi focal stimulation
rng(200)

num_E = '32'; % # of electrodes
file = ['model_cylinder_',num_E,'.mat'];
load(file)

% Define global constraints and maximum current/power values
I_max = 1;
P_max1 = 4;
P_max2 = 5; 

% Initialize weight arrays for targets and non-target regions
w1 = zeros(length(data_x),1);
w2 = zeros(length(data_x),1);
w_non1 = zeros(length(data_x),1);
w_non2 = zeros(length(data_x),1);

% Define target coordinates and their respective radii
target1 = [-8; -6;0]; 
target2 = [0;10;0]; 
r1 = 3.5; 
r2 = 3.5;

% optional, a vector of interest, use when want to know the amplitude in a particular direction
v = [0;1;0];
v = v/norm(v);

% Calculate spatial weights based on distance to target locations
for i = 1:length(w1)
    d1 = xyz(:,i)-target1;
    D1 = norm(d1);
    d2 = xyz(:,i)-target2;
    D2 = norm(d2);
    if norm(xyz(1:2,i))<=25
        if D1 < r1
            w1(i) = (r1-D1+1)^2; 
        else
            w_non1(i) = 1;
        end
        
        if D2 < r2 
            w2(i) = (r2-D2+1)^2; 
        else
            w_non2(i) = 1;
        end
    end
end

% Transfer matrix representing the electric field generated per unit current
A=[data_x;data_y;data_z];
A(isnan(A)) = 0;

% lower and upper bounds
lb =  zeros(4*M,1) ; 
ub = ones(4*M,1) * I_max;

% Initialize random starting points for the optimization solver
a1=rand(M,1); a1= a1*I_max/sum(a1)/M; % f1 
a2=rand(M,1); a2= a2*I_max/sum(a2)/M;
a3=rand(M,1); a3= a3*I_max/sum(a3)/M; % f1 df 
a4=rand(M,1); a4= a4*I_max/sum(a4)/M;
a5=rand(M,1); a5= a5*I_max/sum(a5)/M; % f2  
a6=rand(M,1); a6= a6*I_max/sum(a6)/M;
a7=rand(M,1); a7= a7*I_max/sum(a7)/M; % f2 df 
a8=rand(M,1); a8= a8*I_max/sum(a8)/M;
xi1 = [a1;a2;a3;a4];
xi2 = [a5;a6;a7;a8];

% Define inequality constraints (a*x <= b) for total current limits
a = [ones(1,M) ones(1,M) zeros(1,M) zeros(1,M); zeros(1,M) zeros(1,M) ones(1,M) ones(1,M)];
b = [2*I_max;2*I_max];
% Define equality constraints (Aeq*x = beq) for current conservation (sum = 0)
Aeq = zeros(2,2*M);
Aeq(1,1:M) = 1; Aeq(1,(1+M):2*M) = -1;
Aeq(2,(1+2*M):3*M) = 1; Aeq(2,(1+3*M):4*M) = -1;
beq = [0;0]; 

% Optimization Options
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1000);
options.Display = 'iter-detailed';
options.MaxFunctionEvaluations = 2000000;
options.MaxIterations = 8000;

% Preliminary optimization for Target 1
nonlcon1 =@(x)Pmax_constraint(x,A,w_non1,P_max1,M);
x1 = fmincon(@(X)fun1(X,A,w1,w2,w_non1,M),xi1,a,b,Aeq,beq,lb,ub,nonlcon1,options);
i1_tmp = [x1(1:M)-x1(M+1:2*M);x1(2*M+1:3*M)-x1(3*M+1:4*M)];

% Preliminary optimization for Target 2
nonlcon2 =@(x)Pmax_constraint(x,A,w_non2,P_max1,M);
x2 = fmincon(@(X)fun1(X,A,w2,w1,w_non2,M),xi2,a,b,Aeq,beq,lb,ub,nonlcon2,options); % ,s1,e1
i2_tmp = [x2(1:M)-x2(M+1:2*M);x2(2*M+1:3*M)-x2(3*M+1:4*M)];

% Main Optimization Phase
a = [ones(1,M) ones(1,M) ones(1,M) ones(1,M) zeros(1,M) zeros(1,M) zeros(1,M) zeros(1,M);...
    zeros(1,M) zeros(1,M) zeros(1,M) zeros(1,M) ones(1,M) ones(1,M) ones(1,M) ones(1,M);...
];
b = [4*I_max;4*I_max]; % a*x < b
Aeq = zeros(4,8*M);
Aeq(1,1:M) = 1; Aeq(1,(1+M):2*M) = -1;
Aeq(2,(1+2*M):3*M) = 1; Aeq(2,(1+3*M):4*M) = -1;
Aeq(3,(1+4*M):5*M) = 1; Aeq(3,(1+5*M):6*M) = -1;
Aeq(4,(1+6*M):7*M) = 1; Aeq(4,(1+7*M):8*M) = -1;
beq = [0;0;0;0]; % Aeq*x = beq

nonlcon3 =@(x)Pmax_constraint3(x,A,w_non1&w_non2,P_max2,M,w1,w2);
x3 = fmincon(@(X)fun3(X,A,w1,w2,w_non1&w_non2,M),[x1;x2],a,b,Aeq,beq,[lb;lb],[ub;ub],nonlcon3,options); % ,s1,e1

% Extract optimal current patterns
i_stim1 = [x3(1:M)-x3(M+1:2*M);x3(2*M+1:3*M)-x3(3*M+1:4*M)];
i_stim2 = [x3(1+4*M:M+4*M)-x3(M+1+4*M:2*M+4*M);x3(2*M+1+4*M:3*M+4*M)-x3(3*M+1+4*M:4*M+4*M)];

% Intensity Balancing Phase (Optional)
% Adjust electric field intensity differences between the multiple targets
optional = 0;
if optional
total_I = sum(abs(i_stim1)+abs(i_stim2))/2;
Ex1 = A(1:length(A)/3,:)*i_stim1(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_stim1(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_stim1(1:M);
E1 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_stim1(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_stim1(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_stim1(M+1:2*M);
E2 = [Ex2 Ey2 Ez2];

E_norm1 = vecnorm(E1')';
E_norm2 = vecnorm(E2')';
E_Norm1 = 2*min(abs(E_norm1),abs(E_norm2));
E_max1 = E_MAX(E1,E2);

Ex1 = A(1:length(A)/3,:)*i_stim2(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_stim2(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_stim2(1:M);
E1 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_stim2(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_stim2(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_stim2(M+1:2*M);
E2 = [Ex2 Ey2 Ez2];

E_norm1 = vecnorm(E1')';
E_norm2 = vecnorm(E2')';
E_Norm2 = 2*min(abs(E_norm1),abs(E_norm2));
E_max2 = E_MAX(E1,E2);

if max(E_max1.*(w1>0))<max(E_max2.*(w2>0))
    i_stim1 = max(E_Norm2.*(w2>0))/max(E_Norm1.*(w1>0))*i_stim1;
    i_stim2 = i_stim2;
else
    i_stim2 = max(E_Norm1.*(w1>0))/max(E_Norm2.*(w2>0))*i_stim2;
    i_stim1 = i_stim1;
end
total_I_scale = sum(abs(i_stim1)+abs(i_stim2))/2;
i_stim1 = i_stim1/total_I_scale*2;
i_stim2 = i_stim2/total_I_scale*2;
total_I_scale = sum(abs(i_stim1)+abs(i_stim2))/2;
end
%% Results Visualization
[image1,~] = E_plot_2d(i_stim1,target1,r1,file);
[image2,~] = E_plot_2d(i_stim2,target2,r2,file);
image_max = image1+image2;

figure
axes1 = axes('Parent',gcf);
im = imagesc([-25 25],[-25 25],image_max); % 이거 Y값 부호 반대로 뜨네
colorbar
% title('Emax')
pbaspect([1 1 1])
set(axes1,'YDir','normal');
max_val = max(image_max,[],'all');
clabel = num2str(max_val);
clabel = [clabel(1:4) ' V/m'];
c4= colorbar('Ticks',[0,max_val],'TickLabels',{0,clabel});
hold on
plot(target1(1),target1(2),'r.')
plot(target2(1),target2(2),'r.')
hold on
circle = 0:0.001:2*pi;
plot(25.1*cos(circle),25.1*sin(circle),'w','Linewidth',2)
im.AlphaData = image_max>0;
axis off

%%
function out = fun1(x,A,w1,w2,w_non,M)

i_s = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];

Ex1 = A(1:length(A)/3,:)*i_s(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_s(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_s(1:M);
E1 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_s(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_s(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_s(M+1:2*M);
E2 = [Ex2 Ey2 Ez2];

E_max = E_MAX(E1,E2);
E_max(isnan(E_max)) = 0;

w_other = (w1|w2)-(w1>0);
if sum(w_other) >0
out = -1*w1'*((E_max.*E_max))/sum(w1) +w_other'*(E_max.*E_max)/sum(w_other)/10;
else
out = -w1'*((E_max.*E_max))/sum(w1);
end

end

function out = fun3(x,A,w1,w2,w_non,M)

i_s1 = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];

Ex1 = A(1:length(A)/3,:)*i_s1(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_s1(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_s1(1:M);
E11 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_s1(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_s1(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_s1(M+1:2*M);
E12 = [Ex2 Ey2 Ez2];

E_max1 = E_MAX(E11,E12);
E_max1(isnan(E_max1)) = 0;
% E_norm1 = vecnorm(E1')';
% E_norm2 = vecnorm(E2')';
% E_Norm1 = 2*min(abs(E_norm1),abs(E_norm2));

i_s2 = [x(1+4*M:M+4*M)-x(M+1+4*M:2*M+4*M);x(2*M+1+4*M:3*M+4*M)-x(3*M+1+4*M:4*M+4*M)];

Ex1 = A(1:length(A)/3,:)*i_s2(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_s2(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_s2(1:M);
E21 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_s2(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_s2(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_s2(M+1:2*M);
E22 = [Ex2 Ey2 Ez2];

E_max2 = E_MAX(E21,E22);
E_max2(isnan(E_max2)) = 0;

w_other = (w1<=0)&(w2<=0);
out = -1*(w1'*E_max1/sum(w1)+w2'*E_max2/sum(w2))^2 +w_other'*((E_max1+E_max2).*(E_max1+E_max2))/sum(w_other)/10;

end

function [c,ceq] = Pmax_constraint(x,A,w_non,P_max,M)

i_s = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];
Ex1 = A(1:length(A)/3,:)*i_s(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_s(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_s(1:M);

Ex2 = A(1:length(A)/3,:)*i_s(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_s(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_s(M+1:2*M);

E1 = [Ex1,Ey1,Ez1];
E2 = [Ex2,Ey2,Ez2];
E_max= E_MAX(E1,E2);
E_max(isnan(E_max)) = 0;

c = 1000/sum(w_non)*w_non'*E_max.^2-P_max;
ceq = i_s(1:M).*i_s(M+1:2*M);
end


function [c,ceq] = Pmax_constraint2(x,i_s1,A,w_non,P_max,M)

Ex1 = A(1:length(A)/3,:)*i_s1(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_s1(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_s1(1:M);
E1 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_s1(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_s1(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_s1(M+1:2*M);
E2 = [Ex2 Ey2 Ez2];

E_max1 = E_MAX(E1,E2);
E_max1(isnan(E_max1)) = 0;

i_s2 = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];

Ex1 = A(1:length(A)/3,:)*i_s2(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_s2(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_s2(1:M);
E1 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_s2(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_s2(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_s2(M+1:2*M);
E2 = [Ex2 Ey2 Ez2];

E_max2 = E_MAX(E1,E2);
E_max2(isnan(E_max2)) = 0;


c = 1000/sum(w_non)*w_non'*(E_max1+E_max2).^2-P_max;
ceq = [i_s2(1:M).*i_s2(M+1:2*M)];
end

function [c,ceq] = Pmax_constraint3(x,A,w_non,P_max,M,w1,w2)

i_s1 = [x(1:M)-x(M+1:2*M);x(2*M+1:3*M)-x(3*M+1:4*M)];

Ex1 = A(1:length(A)/3,:)*i_s1(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_s1(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_s1(1:M);
E11 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_s1(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_s1(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_s1(M+1:2*M);
E12 = [Ex2 Ey2 Ez2];

E_max1 = E_MAX(E11,E12);
E_max1(isnan(E_max1)) = 0;

i_s2 = [x(1+4*M:M+4*M)-x(M+1+4*M:2*M+4*M);x(2*M+1+4*M:3*M+4*M)-x(3*M+1+4*M:4*M+4*M)];

Ex1 = A(1:length(A)/3,:)*i_s2(1:M);
Ey1 = A(1+length(A)/3:length(A)*2/3,:)*i_s2(1:M);
Ez1 = A(1+length(A)*2/3:length(A),:)*i_s2(1:M);
E21 = [Ex1 Ey1 Ez1];

Ex2 = A(1:length(A)/3,:)*i_s2(M+1:2*M);
Ey2 = A(1+length(A)/3:length(A)*2/3,:)*i_s2(M+1:2*M);
Ez2 = A(1+length(A)*2/3:length(A),:)*i_s2(M+1:2*M);
E22 = [Ex2 Ey2 Ez2];

E_max2 = E_MAX(E21,E22);
E_max2(isnan(E_max2)) = 0;


w = w1|w2;
P_ = E_max1+E_max2;
c = 1000/sum(w_non)*w_non'*P_.^2-P_max;
ceq =[i_s1(1:M).*i_s1(M+1:2*M); i_s2(1:M).*i_s2(M+1:2*M)];
end