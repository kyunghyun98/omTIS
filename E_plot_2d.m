function [varargout] = E_plot_2d(i_stim,target,r,file)
% function for displaying result (cylinder model)
% i_stim : stimulation current per electrodes and frequency
% target : target point
% r : radius of target area
% file : name of saved data file

load(file,'data_x','data_y','data_z','resolution','M')
n = 25/resolution;
if1=i_stim(1:M); % current parameter for f
if2=i_stim(M+1:2*M); % current parameter for f+df

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

image_x = reshape(envelope_x, [(2*n+1) (2*n+1)])';
image_y = reshape(envelope_y, [(2*n+1) (2*n+1)])';
image_z = reshape(envelope_z, [(2*n+1) (2*n+1)])';
image_max = reshape(E_MAX(E1,E2), [(2*n+1) (2*n+1)])';
image_norm = reshape(envelope_norm, [(2*n+1) (2*n+1)])';
image_x = flip(image_x);
image_y = flip(image_y);
image_z = flip(image_z);
image_max = flip(image_max);
image_norm = flip(image_norm);

max_val = max(E_MAX(E1,E2));
scale = [0,max_val];

point = linspace(0,2*pi,20);
c_x= r*cos(point)+target(1);
c_y= r*sin(point)+target(2);

if 1
    figure
    axes2 = axes('Parent',gcf);
    im = imagesc([-25 25],[-25 25],image_max,scale) ;
    clabel = num2str(max_val);
    clabel = [clabel(1:4) ' V/m'];
    c4= colorbar('Ticks',[0,max_val],'TickLabels',{0,clabel});
    hold on;
    % plot(c_x,c_y)
    pbaspect([1 1 1])
    set(axes2,'YDir','normal');
    title("target = "+num2str(target(1))+", "+num2str(target(2)))
else
    figure
    axes2 = axes('Parent',gcf);
    imagesc([-25 25],[-25 25],image_x,scale)
    colorbar
    title('Ex')
    hold on; plot(c_x,c_y)
    pbaspect([1 1 1])
    set(axes2,'YDir','normal');
    figure
    axes3 = axes('Parent',gcf);
    imagesc([-25 25],[-25 25],image_y,scale)
    colorbar
    title('Ey')
    hold on; plot(c_x,c_y)
    pbaspect([1 1 1])
    set(axes3,'YDir','normal');
    figure
    axes4 = axes('Parent',gcf);
    im = imagesc([-25 25],[-25 25],image_max,scale);
    colorbar
    hold on; plot(c_x,c_y)
    title('Emax')
    pbaspect([1 1 1])
    set(axes4,'YDir','normal');

    figure
    axes1 = axes('Parent',gcf);
    imagesc([-25 25],[-25 25],image_norm,scale)
    colorbar
    hold on; plot(c_x,c_y)
    title('Enorm')
    pbaspect([1 1 1])
    set(axes1,'YDir','normal');
end
im.AlphaData = (flip(image_max)>0);
circle = 0:0.001:2*pi;
plot(25.1*cos(circle),25.1*sin(circle),'w','Linewidth',2)
plot(target(1),target(2),'r.')

varargout{1} = image_max;
varargout{2} = image_norm;


end