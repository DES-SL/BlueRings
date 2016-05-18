function RotPlot(Image_Diff, N_Sq, N_Cir, N_Seg, R_Cir, xcen, ycen, CirMin, CirMax, Center_Type, Im_Cut_Size, PSF_Cut_Size)

FigHandle = figure('Position', [100, 100, 1250, 950], 'Color', [1 1 1]);

if (Center_Type == 1)
    x = linspace((-N_Sq/2)+xcen,(N_Sq/2)+xcen,N_Sq+1);
    y = linspace((-N_Sq/2)+ycen,(N_Sq/2)+ycen,N_Sq+1);
elseif (Center_Type == 2)
    x = linspace((-N_Sq/2)+0.5+xcen,(N_Sq/2)-0.5+xcen,N_Sq);
    y = linspace((-N_Sq/2)+0.5+ycen,(N_Sq/2)-0.5+ycen,N_Sq);
end

imagesc(Image_Diff);
colorbar
hold on
plot(xcen, ycen, 'k.')

ArcHigh = 0.27*(Im_Cut_Size+PSF_Cut_Size);
ArcLow = -0.27*(Im_Cut_Size+PSF_Cut_Size);

set(gca,'XTickMode','manual');
%we will have 3 ticks on X label
set(gca,'XTick',[1, (Im_Cut_Size+PSF_Cut_Size), 2*(Im_Cut_Size+PSF_Cut_Size)]);
set(gca,'XtickLabels',[ArcLow, 0, ArcHigh]);

set(gca,'YTickMode','manual');
%we will have 3 ticks on Y label
set(gca,'YTick',[1, (Im_Cut_Size+PSF_Cut_Size), 2*(Im_Cut_Size+PSF_Cut_Size)]);
set(gca,'YtickLabels',[ArcLow, 0, ArcHigh]);

xlabel('x (arcsec)', 'FontSize', 26); ylabel('y (arcsec)', 'FontSize', 26);
title('Diff Image - Circle Decomposition', 'FontSize', 26)
set(gca,'FontSize',26)

% Horizontal grid 
for k = 1:length(y)
  line([x(1) x(end)], [y(k) y(k)])
end

% Vertical grid
for k = 1:length(y)
  line([x(k) x(k)], [y(1) y(end)])
end

axis square

hold on;

r_min = R_Cir(CirMin);

ang=0:0.01:2*pi; 
xp=r_min*cos(ang);
yp=r_min*sin(ang);
plot(xp+xcen,yp+ycen);

hold on;

r_max = R_Cir(CirMax);

ang=0:0.01:2*pi; 
xp=r_max*cos(ang);
yp=r_max*sin(ang);
plot(xp+xcen,yp+ycen);

for I = 1:N_Seg
    
    ang = ((I-1)/N_Seg)*2*pi;
    
    xp_min = r_min*cos(ang);
    xp_max = r_max*cos(ang);
    
    yp_min = r_min*sin(ang);
    yp_max = r_max*sin(ang);
    
    [xp_min+xcen, xp_max+xcen], [yp_min+ycen, yp_max+ycen]
    
    plot([xp_min+xcen, xp_max+xcen], [yp_min+ycen, yp_max+ycen],'k','LineWidth',2)
    
end