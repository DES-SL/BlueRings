clc;

Image_Blue = fitsread('MockImages/Candidate_1_g_SDSS_img.fits');

%imagesc(Image_Red);

Image_Red = flipud(Image_Red);
Image_Red = fliplr(Image_Red);

[dum, xcen] = max(max(Image_Red));
[dum, ycen] = max(max(Image_Red'));

Image_Blue2 = flipud(Image_Blue);
Image_Blue2 = fliplr(Image_Blue2);
Image_Resi = Image_Blue - Image_Blue2;
imagesc(Image_Resi);

N = 21^2;
x = linspace(-10.5+xcen,10.5+xcen,sqrt(N)+1)
y = linspace(-10.5+ycen,10.5+ycen,sqrt(N)+1)

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

for I = 1:10

r = R_Cir(I+1);

ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(xp+xcen,yp+ycen);

end
r = 1.5;

ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(xp,yp);

r = 2.5;

ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(xp,yp);

r = 3.5;

ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(xp,yp);