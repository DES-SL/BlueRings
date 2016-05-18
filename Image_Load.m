function [Image_i, Image_r, Image_g, Im_xcen, Im_ycen] = Image_Load(Cand, Im_Cut_Size, Plots)

% Load image number 'Cand'

%Image_i = fitsread(['MockImages/Candidate_',num2str(Cand),'_i_SDSS_img.fits']);
%Image_r = fitsread(['MockImages/Candidate_',num2str(Cand),'_r_SDSS_img.fits']);
%Image_g = fitsread(['MockImages/Candidate_',num2str(Cand),'_g_SDSS_img.fits']);

%Image_i = fitsread(['MockImages2/LensingOff_',num2str(Cand),'_i_SDSS_img.fits']);
%Image_r = fitsread(['MockImages2/LensingOff_',num2str(Cand),'_r_SDSS_img.fits']);
%Image_g = fitsread(['MockImages2/LensingOff_',num2str(Cand),'_g_SDSS_img.fits']);
% 
Image_i = fitsread(['MockImages2/LensingOn_',num2str(Cand),'_i_SDSS_img.fits']);
Image_r = fitsread(['MockImages2/LensingOn_',num2str(Cand),'_r_SDSS_img.fits']);
Image_g = fitsread(['MockImages2/LensingOn_',num2str(Cand),'_g_SDSS_img.fits']);

% Image_i = fitsread(['SimPaint/DES2125+0001_3008740963_gal_i.fits']);
% Image_r = fitsread(['SimPaint/DES2125+0001_3008740963_gal_r.fits']);
% Image_g = fitsread(['SimPaint/DES2125+0001_3008740963_gal_g.fits']);

%Center Detection (Improve to look at 3x3 flux squares later

Image_r2 = Image_r;
ImSize = max(size(Image_r2));

CenTot = 10;

xy = zeros(1:CenTot, 1);
Flux = zeros(1:CenTot, 1);
%
for I = 1:CenTot
    %
    [v,ind]=max(Image_r2(:));
    [x,y] = ind2sub(size(Image_r2),ind);
    
    if ((x-1 >= 0) && (y-1 >= 0) && (x+1 <= ImSize) && (y+1 <= ImSize))
        
        Flux(I) = sum(sum(Image_r(x-1:x+1, y-1:y+1)));
        
        Image_r2(x,y) = 0;
        
        xy(I,1) = x;
        xy(I,2) = y;
        
    else
        
        Flux(I) = 0;
    end
    
end
%
[dum, cenin] = max(Flux);
%
Im_xcen = xy(cenin,1);
Im_ycen = xy(cenin,2);

%[dum, Im_xcen] = max(max(Image_r));
%[dum, Im_ycen] = max(max(Image_r'));

% Cut Image given the Cut Size

Im_xmin = Im_xcen-Im_Cut_Size;
Im_xmax = Im_xcen+Im_Cut_Size;
Im_ymin = Im_ycen-Im_Cut_Size;
Im_ymax = Im_ycen+Im_Cut_Size;

%Move center dependent on cuts

Im_xcen = Im_Cut_Size + 1;
Im_ycen = Im_Cut_Size + 1;

%Perform Cut

Image_i = Image_i(Im_xmin:Im_xmax, Im_ymin:Im_ymax);
Image_r = Image_r(Im_xmin:Im_xmax, Im_ymin:Im_ymax);
Image_g = Image_g(Im_xmin:Im_xmax, Im_ymin:Im_ymax);

% if strcmp(Plots, 'On')

%ArcLow = -0.27*(Im_Cut_Size/2);
%ArcHigh = 0.27*(Im_Cut_Size/2);

% 
% % First Plot the 'red' Image
% 
% subplot(3,3,1);
% imagesc(Image_i);
% hold on
% colorbar
% plot(Im_xcen, Im_ycen, 'k.')
% 
% set(gca,'XTickMode','manual');
% %we will have 3 ticks on X label
% set(gca,'XTick',[1, Im_Cut_Size, 2*Im_Cut_Size]);
% set(gca,'XtickLabels',[ArcLow, 0, ArcHigh]);
% 
% set(gca,'YTickMode','manual');
% %we will have 3 ticks on Y label
% set(gca,'YTick',[1, Im_Cut_Size, 2*Im_Cut_Size]);
% set(gca,'YtickLabels',[ArcLow, 0, ArcHigh]);
% 
% xlabel('x (arcsec)', 'FontSize', 15); ylabel('y (arcsec)', 'FontSize', 15);
% title('Image 1 (Redder Band)', 'FontSize', 15)
% set(gca,'FontSize',16)
% 
% 
% % Next plot the 'Blue' Image
% 
% subplot(3,3,2);
% imagesc(Image_g);
% hold on
% colorbar
% plot(Im_xcen, Im_ycen, 'k.')
% 
% set(gca,'XTickMode','manual');
% %we will have 3 ticks on X label
% set(gca,'XTick',[1, Im_Cut_Size, 2*Im_Cut_Size]);
% set(gca,'XtickLabels',[ArcLow, 0, ArcHigh]);
% 
% set(gca,'YTickMode','manual');
% %we will have 3 ticks on Y label
% set(gca,'YTick',[1, Im_Cut_Size, 2*Im_Cut_Size]);
% set(gca,'YtickLabels',[ArcLow, 0, ArcHigh]);
% 
% xlabel('x (arcsec)', 'FontSize', 15); ylabel('y (arcsec)', 'FontSize', 15);
% title('Image 2 (Blue Band)', 'FontSize', 15)
% set(gca,'FontSize',16)
% 
% end
% 













