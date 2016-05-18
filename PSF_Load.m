function [PSF_i, PSF_r, PSF_g, PSF_xcen, PSF_ycen] = PSF_Load(Cand, PSF_Cut_Size)

% Load PSF and Cut

%PSF_i = fitsread(['MockImages/Candidate_',num2str(Cand),'_i_SDSS_psf.fits']);
%PSF_r = fitsread(['MockImages/Candidate_',num2str(Cand),'_r_SDSS_psf.fits']);
%PSF_g = fitsread(['MockImages/Candidate_',num2str(Cand),'_g_SDSS_psf.fits']);

%PSF_i = fitsread(['MockImages2/LensingOff_',num2str(Cand),'_i_SDSS_psf.fits']);
%PSF_r = fitsread(['MockImages2/LensingOff_',num2str(Cand),'_r_SDSS_psf.fits']);
%PSF_g = fitsread(['MockImages2/LensingOff_',num2str(Cand),'_g_SDSS_psf.fits']);
% 
PSF_i = fitsread(['MockImages2/LensingOn_',num2str(Cand),'_i_SDSS_psf.fits']);
PSF_r = fitsread(['MockImages2/LensingOn_',num2str(Cand),'_r_SDSS_psf.fits']);
PSF_g = fitsread(['MockImages2/LensingOn_',num2str(Cand),'_g_SDSS_psf.fits']);

% PSF_i = fitsread(['SimPaint/DES2125+0001_3008740963_psf_i.fits']);
% PSF_r = fitsread(['SimPaint/DES2125+0001_3008740963_psf_r.fits']);
% PSF_g = fitsread(['SimPaint/DES2125+0001_3008740963_psf_g.fits']);

% Cut PSF

%NEED IF LOOP MAKING CUT SIZE .5 IF Center is a half integer, FURTHER Need
%Cutsize X, Y for caseswhere only one is .5

PSF_xcen = 38.5;
PSF_ycen = 38.5;

PSF_xmin = PSF_xcen-PSF_Cut_Size;
PSF_xmax = PSF_xcen+PSF_Cut_Size;
PSF_ymin = PSF_ycen-PSF_Cut_Size;
PSF_ymax = PSF_ycen+PSF_Cut_Size;

PSF_xcen = PSF_Cut_Size + 1;
PSF_ycen = PSF_Cut_Size + 1;

PSF_i = PSF_i(PSF_xmin:PSF_xmax, PSF_ymin:PSF_ymax);
PSF_r = PSF_r(PSF_xmin:PSF_xmax, PSF_ymin:PSF_ymax);
PSF_g = PSF_g(PSF_xmin:PSF_xmax, PSF_ymin:PSF_ymax);
