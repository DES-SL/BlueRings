%%Rot Code


Sub_Type = 'Rot';

Image_g2 = imrotate(Image_g, 180);
Image_Diff = Image_g% - Image_g2;

%[dum, xcen] = max(max(Image_g2));
%[dum, ycen] = max(max(Image_g2'));

if ((ceil(xcen) - xcen) == 0) && ((ceil(ycen) - ycen) == 0)
    Center_Type = 1;
elseif ((ceil(xcen) - xcen) == 0.5) && ((ceil(ycen) - ycen) == 0.5)
    Center_Type = 2;
else
    Center_Type = 3;
end

xcen = Im_xcen;
ycen = Im_ycen;

%DiffPlot(Image_i, Image_g, Image_Diff, N_Sq, N_Cir, R_Cir, xcen, ycen, Center_Type)

%imagesc(Image_Diff)
