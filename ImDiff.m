function [Image_Diff, xcen, ycen, Center_Type, ImDiffMin] = ImDiff(Image_1, Image_2, PSF_1, PSF_2, PSF_xcen, PSF_ycen, Im_xcen, Im_ycen, R_Min, R_Max, Plots)

%Convert center of image (Found with flux summing) to center of PSF
%convolved Images

xcen = PSF_xcen + Im_xcen - 1.0;
ycen = PSF_ycen + Im_ycen - 1.0;

%Concolve image 1 with 2, visa versa

Im_1_conv_2 = conv2(Image_1, PSF_2);
Im_2_conv_1 = conv2(Image_2, PSF_1);

% Alpha MAsk - Only pixels in this mask are used to calc alpha
Alpha_Mask = zeros(max(size(Im_1_conv_2)), max(size(Im_1_conv_2)));

for I = 1:max(size(Im_1_conv_2))
    for J = 1:max(size(Im_2_conv_1))
        x_arc = 0.27*(I-xcen);
        y_arc = 0.27*(J-ycen);
        r_arc = sqrt(x_arc^2 + y_arc^2);
        if (r_arc > R_Min && r_arc < R_Max)
            Alpha_Mask(I,J) = 1;
        end
    end
end

alphasum = 0;
tot = 0;

for I = 1:max(size(Im_1_conv_2))
    for J = 1:max(size(Im_2_conv_1))
        if (Alpha_Mask(I,J) == 1)
        alphasum = alphasum + (Im_1_conv_2(I,J) / Im_2_conv_1(I,J));
        tot = tot+1;
        end
    end
end

alpha = alphasum/tot;

% perform difference image, subtracting image 1 from 2*alpha (both
% convolved)

Image_Diff = Im_1_conv_2(:,:) - alpha*Im_2_conv_1(:,:);

%Center Type:
%Is the center in the middle of a pixel (=1)
% Between four pixels on a corner (=2)
% Between two pixels (=3)

if ((ceil(xcen) - xcen) == 0) && ((ceil(ycen) - ycen) == 0)
    Center_Type = 1;
elseif ((ceil(xcen) - xcen) == 0.5) && ((ceil(ycen) - ycen) == 0.5)
    Center_Type = 2;
else
    Center_Type = 3;
end

ImDiffMin = min(min(Image_Diff));

% 
%ImDiffCut = 19;

%Image_Diff(:, 1:ycen-ImDiffCut) = 0;
%Image_Diff(:, ycen+ImDiffCut:60) = 0;
%Image_Diff(1:xcen-ImDiffCut, :) = 0;
%Image_Diff(xcen+ImDiffCut:60,:) = 0;

% FigHandle = figure('Position', [100, 100, 1250, 950], 'Color', [1 1 1]);

%Now plot the image difference

% imagesc(Im_1_conv_2);
% colorbar
% hold on
% 
% Im_Cut_Size = 21;
%     
%     PSF_Cut_Size = 9.5;
% 
% ArcHigh = 0.27*(Im_Cut_Size+PSF_Cut_Size);
% ArcLow = -0.27*(Im_Cut_Size+PSF_Cut_Size);
% 
% set(gca,'XTickMode','manual');
% %we will have 3 ticks on X label
% set(gca,'XTick',[1, (Im_Cut_Size+PSF_Cut_Size), 2*(Im_Cut_Size+PSF_Cut_Size)]);
% set(gca,'XtickLabels',[ArcLow, 0, ArcHigh]);
% 
% set(gca,'YTickMode','manual');
% %we will have 3 ticks on Y label
% set(gca,'YTick',[1, (Im_Cut_Size+PSF_Cut_Size), 2*(Im_Cut_Size+PSF_Cut_Size)]);
% set(gca,'YtickLabels',[ArcLow, 0, ArcHigh]);
% 
% xlabel('x (arcsec)', 'FontSize', 26); ylabel('y (arcsec)', 'FontSize', 26);
% title('Image I Conv w/ PSF G', 'FontSize', 26)
% set(gca,'FontSize',26)

% if strcmp(Plots, 'On')
%     
% ArcHigh = ImDiffCut*0.27;
% ArcLow = -ImDiffCut*0.27;
% 
% subplot(3,3,3)
% 
% Image_Diff_Plot = Image_Diff((ceil(xcen)-ImDiffCut):(floor(xcen)+ImDiffCut), (ceil(ycen)-ImDiffCut):(floor(ycen)+ImDiffCut));
% 
% %Now plot the image difference
% 
% imagesc(Image_Diff_Plot);
% colorbar
% hold on
% plot(xcen, ycen, 'k.')
% 
% set(gca,'XTickMode','manual');
% %we will have 3 ticks on X label
% set(gca,'XTick',[1, ImDiffCut, 2*ImDiffCut]);
% set(gca,'XtickLabels',[ArcLow, 0, ArcHigh]);
% 
% set(gca,'YTickMode','manual');
% %we will have 3 ticks on Y label
% set(gca,'YTick',[1, ImDiffCut, 2*ImDiffCut]);
% set(gca,'YtickLabels',[ArcLow, 0, ArcHigh]);
% 
% xlabel('x (arcsec)', 'FontSize', 15); ylabel('y (arcsec)', 'FontSize', 15);
% title('Difference Image', 'FontSize', 15)
% set(gca,'FontSize',16)
% 
% end
end
