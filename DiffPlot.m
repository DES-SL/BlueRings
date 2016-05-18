function DiffPlot(Image_i, Image_r, Image_g, Image_Diff, N_Sq, N_Cir, R_Cir, xcen, ycen, Center_Type, Im_Cut_Size, PSF_Cut_Size)

if (Center_Type == 1)
    x = linspace((-N_Sq/2)+xcen,(N_Sq/2)+xcen,N_Sq+1);
    y = linspace((-N_Sq/2)+ycen,(N_Sq/2)+ycen,N_Sq+1);
elseif (Center_Type == 2)
    x = linspace((-N_Sq/2)+0.5+xcen,(N_Sq/2)-0.5+xcen,N_Sq);
    y = linspace((-N_Sq/2)+0.5+ycen,(N_Sq/2)-0.5+ycen,N_Sq);
end

FigHandle = figure('Position', [100, 100, 1250, 950], 'Color', [1 1 1]);

%Now plot the image difference

imagesc(Image_r);
colorbar
hold on
plot(xcen, ycen, 'k.')

ArcHigh = 0.27*(Im_Cut_Size);
ArcLow = -0.27*(Im_Cut_Size);

set(gca,'XTickMode','manual');
%we will have 3 ticks on X label
set(gca,'XTick',[1, (Im_Cut_Size), 2*(Im_Cut_Size)]);
set(gca,'XtickLabels',[ArcLow, 0, ArcHigh]);

set(gca,'YTickMode','manual');
%we will have 3 ticks on Y label
set(gca,'YTick',[1, (Im_Cut_Size), 2*(Im_Cut_Size)]);
set(gca,'YtickLabels',[ArcLow, 0, ArcHigh]);

xlabel('x (arcsec)', 'FontSize', 26); ylabel('y (arcsec)', 'FontSize', 26);
title('Image (R band)', 'FontSize', 26)
set(gca,'FontSize',26)

FigHandle = figure('Position', [100, 100, 1250, 950], 'Color', [1 1 1]);

%Now plot the image difference

imagesc(Image_i);
colorbar
hold on
plot(xcen, ycen, 'k.')

ArcHigh = 0.27*(Im_Cut_Size);
ArcLow = -0.27*(Im_Cut_Size);

set(gca,'XTickMode','manual');
%we will have 3 ticks on X label
set(gca,'XTick',[1, (Im_Cut_Size), 2*(Im_Cut_Size)]);
set(gca,'XtickLabels',[ArcLow, 0, ArcHigh]);

set(gca,'YTickMode','manual');
%we will have 3 ticks on Y label
set(gca,'YTick',[1, (Im_Cut_Size), 2*(Im_Cut_Size)]);
set(gca,'YtickLabels',[ArcLow, 0, ArcHigh]);

xlabel('x (arcsec)', 'FontSize', 26); ylabel('y (arcsec)', 'FontSize', 26);
title('Image (I band)', 'FontSize', 26)
set(gca,'FontSize',26)

FigHandle = figure('Position', [100, 100, 1250, 950], 'Color', [1 1 1]);

%Now plot the image difference

imagesc(Image_g);
colorbar
hold on
plot(xcen, ycen, 'k.')

ArcHigh = 0.27*(Im_Cut_Size);
ArcLow = -0.27*(Im_Cut_Size);

set(gca,'XTickMode','manual');
%we will have 3 ticks on X label
set(gca,'XTick',[1, (Im_Cut_Size), 2*(Im_Cut_Size)]);
set(gca,'XtickLabels',[ArcLow, 0, ArcHigh]);

set(gca,'YTickMode','manual');
%we will have 3 ticks on Y label
set(gca,'YTick',[1, (Im_Cut_Size), 2*(Im_Cut_Size)]);
set(gca,'YtickLabels',[ArcLow, 0, ArcHigh]);

xlabel('x (arcsec)', 'FontSize', 26); ylabel('y (arcsec)', 'FontSize', 26);
title('Image (G band)', 'FontSize', 26)
set(gca,'FontSize',26)


FigHandle = figure('Position', [100, 100, 1250, 950], 'Color', [1 1 1]);

%Now plot the image difference

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
title('Difference Image (-ve = Blue Residuals)', 'FontSize', 26)
set(gca,'FontSize',26)

% % Horizontal grid 
% for k = 1:length(y)
%   line([x(1) x(end)], [y(k) y(k)])
% end
% 
% % Vertical grid
% for k = 1:length(y)
%   line([x(k) x(k)], [y(1) y(end)])
% end
% 
% axis square
% 
% hold on;
% 
% for I = 1:N_Cir
% 
% r = R_Cir(I);
% 
% ang=0:0.01:2*pi; 
% xp=r*cos(ang);
% yp=r*sin(ang);
% plot(xp+xcen,yp+ycen);
% 
% end

end