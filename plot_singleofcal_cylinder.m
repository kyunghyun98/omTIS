% Ploting simulated result of stimulation parameters for single focal
% stimulation in shpere model
% Manual setting or loading data for stimulation current
mannual = 0;
if mannual ==1
    M = 16;
    i_f1 = zeros(M,1);
    i_f2 = zeros(M,1);

    i_f1(1) =1;
    i_f1(2) =-1;
    i_f2(8) =1;
    i_f2(9) =-1;
    i_stim = [i_f1; i_f2];
else
    load("Fig2\Fig2c_1.mat","i_stim","M")
end
%%
file = ['model_cylinder_',num2str(M),'.mat'];
load(file,'data_x','data_y','data_z','resolution','M')
n = 25/resolution;
if1=i_stim(1:M);
if2=i_stim(M+1:2*M);

E1x = data_x*if1;
E1y = data_y*if1;
E1z = data_z*if1;

E2x = data_x*if2;
E2y = data_y*if2;
E2z = data_z*if2;

E1 = [E1x,E1y,E1z];
E2 = [E2x,E2y,E2z];

envelope_x = 2*min(abs(E1x),abs(E2x));
envelope_y = 2*min(abs(E1y),abs(E2y));
envelope_z = 2*min(abs(E1z),abs(E2z));

E_norm1 = vecnorm(E1')';
E_norm2 = vecnorm(E2')';
envelope_norm = 2*min(abs(E_norm1),abs(E_norm2));

image_x = flip(reshape(envelope_x, [(2*n+1) (2*n+1)])');
image_y = flip(reshape(envelope_y, [(2*n+1) (2*n+1)])');
image_z = flip(reshape(envelope_z, [(2*n+1) (2*n+1)])');
image_max = flip(reshape(E_MAX(E1,E2), [(2*n+1) (2*n+1)])');
image_norm = flip(reshape(envelope_norm, [(2*n+1) (2*n+1)])');

max_val = max(E_MAX(E1,E2));
scale = [0,max_val];

figure
axes2 = axes('Parent',gcf);
im = imagesc([-25 25],[-25 25],image_max,scale) ;
clabel = num2str(max_val);
clabel = [clabel ' V/m'];
c4= colorbar('Ticks',[0,max_val],'TickLabels',{0,clabel});
hold on;
pbaspect([1 1 1])
set(axes2,'YDir','normal');

im.AlphaData = (image_max)>0;
circle = 0:0.001:2*pi;
plot(25.1*cos(circle),25.1*sin(circle),'w','Linewidth',2)

