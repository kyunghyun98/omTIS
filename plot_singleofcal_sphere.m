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
    load("Fig2\Fig2f.mat","i_stim","M")
end
%%
file = ['model_sphere_',num2str(M),'.mat'];
load(file,'data_x','data_y','data_z','resolution','M','R','xyz')
n = R/resolution;
x_num = (target(1)+R)/resolution + 1;
y_num = (target(2)+R)/resolution + 1;
z_num = (target(3)+R)/resolution + 1;

if1=i_stim(1:M);
if2=i_stim(M+1:M*2);

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

image3d_x = reshape(envelope_x, [(2*n+1) (2*n+1) (2*n+1)]);
image3d_y = reshape(envelope_y, [(2*n+1) (2*n+1) (2*n+1)]);
image3d_z = reshape(envelope_z, [(2*n+1) (2*n+1) (2*n+1)]);

image3d_max = reshape(E_MAX(E1,E2), [(2*n+1) (2*n+1) (2*n+1)]);
distance = reshape(xyz.distance, [(2*n+1) (2*n+1) (2*n+1)]);

image_tmp = image3d_max;
image2d{1} = reshape(image_tmp(z_num,:,:),[(2*n+1),(2*n+1)]);
image2d{2} = reshape(image_tmp(:,y_num,:),[(2*n+1),(2*n+1)]);
image2d{3} = reshape(image_tmp(:,:,x_num),[(2*n+1),(2*n+1)]);

distance_2d{1} = reshape(distance(z_num,:,:),[(2*n+1),(2*n+1)]);
distance_2d{2} = reshape(distance(:,y_num,:),[(2*n+1),(2*n+1)]);
distance_2d{3} = reshape(distance(:,:,x_num),[(2*n+1),(2*n+1)]);
scale = [0 max(image_tmp.*(distance<=25*0.842),[],'all')];

tiledlayout(2,2);
ax1 = nexttile;
[X,Y,Z] = sphere(ax1);
title('model')
X = X * R;
Y = Y * R;
Z = Z * R;
surf(X,Y,Z,zeros(size(X)),'FaceAlpha',0.1,'FaceColor',[0.3010 0.7450 0.9330] ,'EdgeColor','none')
map = [0.8 0.8 0.8];

hold on
[x,y] = meshgrid(-R:0.5:R,-R:0.5:R);
planecolor = [0.5 0.5 0.5];
for i = 1:3
    switch i
        case 1
            sf =surf(x,y,zeros(size(x))+target(4-i),image2d{i}.*(distance_2d{i}<=25*0.842),'FaceAlpha','flat','EdgeColor','none');
            alpha_map = (distance_2d{i}<=25*0.842);
            sf.AlphaData = alpha_map*0.8;
        case 2
            sf =surf(x,zeros(size(x))+target(4-i),y,image2d{i}.*(distance_2d{i}<=25*0.842) ,'FaceAlpha','flat','EdgeColor','none');
            alpha_map = (distance_2d{i}<=25*0.842);
            sf.AlphaData = alpha_map*0.8;
        case 3
            sf = surf(zeros(size(x))+target(4-i),x,y,image2d{i}.*(distance_2d{i}<=25*0.842),'FaceAlpha','flat','EdgeColor','none');
            alpha_map = (distance_2d{i}<=25*0.842);
            sf.AlphaData = alpha_map*0.8;
    end
end
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title("target = ("+num2str(target(1))+", "+num2str(target(2))+", "+num2str(target(3))+")")
axis equal
hold off

point = linspace(0,2*pi,20);
idx1 = ['zyx';'yzx';'xzy'];
idx2 = [1 2;1 3;2 3];

for i = 1:3
    axes1 = nexttile;
    max_val = max(image2d{i}.*(distance_2d{i}<=25*0.842),[],'all');

    im = imagesc([-R R],[-R R],image2d{i}.*(distance_2d{i}<=25*0.842),scale);
    alpha_map = (distance_2d{i}<=25*0.842);
    im.AlphaData = alpha_map;
    c1 = colorbar;

    max_val = scale(2);
    clabel = num2str(round(max_val*1000)/1000);
    clabel = [clabel ' V/m'];
    c1= colorbar('Ticks',[0,max_val],'TickLabels',{0,clabel});
    hold on

    % plot(target(1),target(2),'r.')
    title(['Plane ', idx1(i,1),' = ',num2str(target(4-i))])

    pbaspect([1 1 1])
    ylabel([idx1(i,2) ' (mm)'])
    xlabel([idx1(i,3) ' (mm)'])
    set(axes1,'YDir','normal');
end