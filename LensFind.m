%clear all; clc;

%% MC init

% First Define set of circles

N_Cir = 16;

Xc_Cir = 0;
Yc_Cir = 0;

R_Cir_Space = 1.5;

% Number of Circle Segments

N_Seg = 32;

%Now Define Square Grid

N_Sq = 41;

% Number of MC iterations?

N_Tot = 75000;

Plots = 'Off';
Plots = 'On';

%Comment this out if you have the areas as can be time-consuming to compute
[Area_Frac, Area_Cir, Cir_Cross, R_Cir, R_Pix_Max, R_Pix_Min, Seg_Frac, Seg_Grid, Seg_Theta] = MC_Areas(N_Cir, R_Cir_Space, Xc_Cir, Yc_Cir, N_Sq, N_Seg, N_Tot);

%% Load images, PSF, centers and Cut

Cand_Index = 0;

tic

% Candidates used for Analysis
for Cand = 1:1
   
    
    %% Image / PSF Loading
    
    % The Cut Size trims the image around te center by cut size each side and
    % then loads the image
    
    Im_Cut_Size = 21;
    [Image_i, Image_r, Image_g, Im_xcen, Im_ycen] = Image_Load(Cand, Im_Cut_Size, Plots);
    
    % The Cut Size trims the image around te center by cut size each side and
    % then loads the PSF
    
    PSF_Cut_Size = 9.5;
    [PSF_i, PSF_r, PSF_g, PSF_xcen, PSF_ycen] = PSF_Load(Cand, PSF_Cut_Size);
    
    
    %% Difference Image / Rotate Subtraction
    
    % Cuts which Circle (Res Cut 1) and Segment (Res Cut 2) must be above
    
    Res_Cut1 = -1.2;
    Res_Cut2 = -3.;
    
    % Radii which alpha, image diff scale, is calculated on
    
    R_Min = 0.5;
    R_Max = 2.0;
    
    [Image_Diff, xcen, ycen, Center_Type, ImDiffMin] = ImDiff(Image_i, Image_g, PSF_i, PSF_g, PSF_xcen, PSF_ycen, Im_xcen, Im_ycen, R_Min, R_Max, Plots);
    %[Image_Diff, xcen, ycen, Center_Type] = ImDiff(Image_r, Image_g, PSF_r, PSF_g, PSF_xcen, PSF_ycen, Im_xcen, Im_ycen, R_Min, R_Max, Plots);
    
    Res_Min1 = min([ImDiffMin*0.1, Res_Cut1]);
    Res_Min2 = min([ImDiffMin*0.2, Res_Cut2]);
    
    % PLot Image DIff with grid / circles
    if strcmp(Plots, 'On')
        DiffPlot(Image_i, Image_r, Image_g, Image_Diff, N_Sq, N_Cir, R_Cir, xcen, ycen, Center_Type, Im_Cut_Size, PSF_Cut_Size);
    end
    
    
    %% Calc Flux in concentric circles
    
    CircFlux = zeros(1,N_Cir);
    
    % Depending on the center, x / y start and end will go from slightly
    % different pixels
    
    if (Center_Type == 1)
        
        x_start = xcen - N_Sq/2 + 0.5;
        x_end = xcen + N_Sq/2 - 0.5;
        y_start = ycen - N_Sq/2 + 0.5;
        y_end = ycen + N_Sq/2 - 0.5;
        
    elseif (Center_Type == 2)
        
        x_start = xcen - N_Sq/2 + 1;
        x_end = xcen + N_Sq/2 - 1;
        y_start = ycen - N_Sq/2 + 1;
        y_end = ycen + N_Sq/2 - 1;
        
    end
    
    % Go to each square pixel, check it has negative (Blue) residual flux
    % If so, use the Area fractiosn to put that flux into the appropriate
    % circle.
    
    CircFlux = zeros(N_Cir,1);
    
    for I = x_start:x_end
        for J = y_start:y_end
            if (Image_Diff(I,J) < 0)
                for K = 1:N_Cir
                    CircFlux(K) = CircFlux(K) + (Image_Diff(I,J))*Area_Frac(I-x_start+1,J-y_start+1,K, Center_Type);
                end
            end
        end
    end
    
    % Pick the circles which contain significant blue residuals
    
    clear CircUse
    In = 1;
    Flux_Detect = 0;
    
    for K = 1:N_Cir
        CircFlux(K) = CircFlux(K)/Area_Cir(K);
        if (CircFlux(K) < Res_Min1)
            
            % Inner flux is often ignored due to small circle... if circle 3 makes the
            % cut then also include circle 2
            
            if K == 3
                
                CircUse(In) = K-1;
                In = In + 1;
            end
            
            CircUse(In) = K;
            In = In + 1;
            Flux_Detect = 1;
        end
    end
    
    %% Segment Area Splitting
    
    % Flux Detect is switched to 1 if blue residuals in the circles are found -
    % if no blue residuals then no point splitting in to segments
    
    if (Flux_Detect == 1)
        
        % Biggest / Smallest Circles and their radii
        
        CirMin = min(CircUse);
        CirMax = max(CircUse);
        
        R_CirMin = R_Cir(CirMin);
        R_CirMax = R_Cir(CirMax);
        
        %Plot the Circles with segment cuts
        
        if strcmp(Plots, 'On')
            
            RotPlot(Image_Diff, N_Sq, N_Cir, N_Seg, R_Cir, xcen, ycen, CirMin, CirMax, Center_Type, Im_Cut_Size, PSF_Cut_Size)
            
        end
        
        % Now calc flux in all segments
        
        Seg_Flux = zeros(N_Seg,1);
        Seg_Flux2 = zeros(N_Seg,1);
        
        for I = x_start:x_end
            for J = y_start:y_end
                
                if (Image_Diff(I,J) < 0)
                    
                    x_Pix = I-x_start-(N_Sq/2);
                    y_Pix = J-y_start-(N_Sq/2);
                    
                    % Is square pixel within the circular annuli?
                    
                    if (R_Pix_Max(I-x_start+1, J-y_start+1, Center_Type) > R_CirMin && R_Pix_Min(I-x_start+1, J-y_start+1, Center_Type) < R_CirMax);
                        
                        %Now decide if pixel is on an inside circle edge - if
                        %it is draw the correct Cir_fraction value
                        
                        if (R_Pix_Min(I-x_start+1, J-y_start+1, Center_Type) < R_CirMin && R_Pix_Max(I-x_start+1, J-y_start+1, Center_Type) > R_CirMin);
                            
                            Cir_Frac = 1 - Area_Frac(I-x_start+1, J-y_start+1, CirMin, Center_Type);
                            
                            % and an Outside Circle edge - again getting the
                            % appropriate fraction if so
                            
                        elseif (R_Pix_Min(I-x_start+1, J-y_start+1, Center_Type) < R_CirMax && R_Pix_Max(I-x_start+1, J-y_start+1, Center_Type) > R_CirMax);
                            
                            Cir_Frac = Area_Frac(I-x_start+1, J-y_start+1, CirMax, Center_Type) + Area_Frac(I-x_start+1, J-y_start+1, CirMax-1, Center_Type);
                            
                        else
                            
                            % Else the pixel is not on a circle edge and is
                            % fully inside a circle
                            
                            Cir_Frac = 1;
                            
                        end
                        
                        % The Segnment number of each square pixels were
                        % precalculated on Seg_Grid
                        
                        Seg_No =  Seg_Grid(I-x_start+1, J-y_start+1, Center_Type);
                        
                        % If Seg Grid has a decimal point (e.g. Seg_No = 1.5)
                        % then that means the square pixel is split accross two
                        % segments and appropriate area fractions must be used
                        % to split the square pixel flux
                        
                        if ( mod(Seg_No,1) == 0 )
                            
                            % No segment split, use full area
                            
                            Seg_Flux(Seg_No) = Seg_Flux(Seg_No) + Image_Diff(I,J)*Cir_Frac;
                            
                        elseif ( mod(Seg_No,1) == 0.5)
                            
                            % Segmnet splits square pixel, split flux given the
                            % areas calculated in Seg_Frac
                            
                            if (floor(Seg_No) == 0)
                                Seg_No = N_Seg;
                            end
                            
                            Seg_Flux(floor(Seg_No)) = Seg_Flux(floor(Seg_No)) + Image_Diff(I,J)*Cir_Frac*Seg_Frac(I-x_start+1,J-y_start+1, Center_Type);
                            
                            Seg_Flux(ceil(Seg_No)) = Seg_Flux(ceil(Seg_No)) + Image_Diff(I,J)*Cir_Frac*(1 - Seg_Frac(I-x_start+1,J-y_start+1, Center_Type));
                            
                        end
                        
                    end
                end
            end
        end
        
        % Calculate Segment areas
        Seg_Area = ((pi*(R_CirMax^2)) - (pi*(R_CirMin^2))) / N_Seg;
        
        % Reset image subtend / flux calcs and counting vairables
        Seg_ImSubtend(1,1:N_Seg) = 0;
        Seg_ImFlux(1,1:N_Seg) = 0;
        ImNo = 1;
        Counting = 0;
        
        %Divide fluxes by area
        
        for I = 1:N_Seg
            Seg_Flux(I) = Seg_Flux(I) / Seg_Area;
        end
        
        K = 1;
        
        %First reorder Seg_Flux, so that it begins rotating at the highest
        %residual value (So that images arn't split due to a poor choice of
        %starting point
        
        [dum, index] = max(Seg_Flux);
        
        %         for I = index:N_Seg
        %             Seg_Flux2(K) = Seg_Flux(I);
        %             K = K+1;
        %         end
        %
        %         if index > 1
        %             for I = 1:index-1
        %                 Seg_Flux2(K) = Seg_Flux(I);
        %                 K = K+1;
        %             end
        %         end
        
        % With the flux reordered, now scan the segments, picking out those
        % below a minimum flux cut. When this happens, begin 'counting' a
        % continuous image.
        
        % For each image this gives the ImSubtend - the angle subtended by the
        % blue residual.
        
        % If an image then ends (part of the image subtend falls below the
        % cut), a second image will be detected if a new segement goes below
        % the residual cut
        
        Seg_Plots = zeros(1,N_Seg);
        
        clear Seg_ImAngStartEnd Seg_ImAngAvg
        
        K = 0;
        
        for I = index:N_Seg
            
            K = K+1;
            
            Seg_Flux2(K) = Seg_Flux(I);
            Seg_Theta2(K) = K*(360/N_Seg);
            
            if (Seg_Flux(I) < Res_Min2)
                
                % Store starting angle
                
                if Counting == 0
                    Seg_ImAngStartEnd(ImNo,1) = K*(360/(N_Seg));
                end
                
                
                % Counting 2 is used to merge two segments, where one is under
                % the cut but the next is then above. Thus this if loops
                % includes that missed segements flux, angle.
                
                if Counting == 2
                    
                    Seg_ImSubtend(ImNo) = Seg_ImSubtend(ImNo) + (360/N_Seg);
                    Seg_ImFlux(ImNo) = Seg_ImFlux(ImNo) + Seg_Flux(I-1);
                    Seg_Plots(I-1) = 1;
                    
                end
                
                Counting = 1;
                
                Seg_ImSubtend(ImNo) = Seg_ImSubtend(ImNo) + (360/N_Seg);
                Seg_ImFlux(ImNo) = Seg_ImFlux(ImNo) + Seg_Flux(I);
                
                Seg_Plots(I) = 1;
                
                
            elseif (Seg_Flux(I) > Res_Min2) && Counting == 1
                
                Counting = 2;
                
            elseif (Seg_Flux(I) > Res_Min2) && Counting == 2
                
                Seg_ImAngStartEnd(ImNo,2) = (K-1)*(360/N_Seg);
                
                Counting = 0;
                Seg_ImAng(1) = 0;
                
                if ImNo > 1
                    Seg_ImAng(ImNo) = Seg_ImSubtend(ImNo);
                end
                
                ImNo = ImNo + 1;
                
            end
            
        end
        
        for I = 1:index-1
            
                 K = K+1;
                 
            Seg_Flux2(K) = Seg_Flux(I);
            Seg_Theta2(K) = K*(360/N_Seg);
            
            
            if (Seg_Flux(I) < Res_Min2)
                
                % Store starting angle
                
                if Counting == 0
                    Seg_ImAngStartEnd(ImNo,1) = K*(360/(N_Seg));
                end
                
                
                % Counting 2 is used to merge two segments, where one is under
                % the cut but the next is then above. Thus this if loops
                % includes that missed segements flux, angle.
                
                if Counting == 2
                    
                    Seg_ImSubtend(ImNo) = Seg_ImSubtend(ImNo) + (360/N_Seg);
                    Seg_ImFlux(ImNo) = Seg_ImFlux(ImNo) + Seg_Flux(I-1);
                    Seg_Plots(I-1) = 1;
                    
                end
                
                Counting = 1;
                
                Seg_ImSubtend(ImNo) = Seg_ImSubtend(ImNo) + (360/N_Seg);
                Seg_ImFlux(ImNo) = Seg_ImFlux(ImNo) + Seg_Flux(I);
                
                Seg_Plots(I) = 1;
                
                
            elseif (Seg_Flux(I) > Res_Min2) && Counting == 1
                
                Counting = 2;
                
            elseif (Seg_Flux(I) > Res_Min2) && Counting == 2
                
                Seg_ImAngStartEnd(ImNo,2) = (K-1)*(360/N_Seg);
                
                Counting = 0;
                Seg_ImAng(1) = 0;
                
                if ImNo > 1
                    Seg_ImAng(ImNo) = Seg_ImSubtend(ImNo);
                end
                
                ImNo = ImNo + 1;
                
            end
            
        end
        
              
        ImNo = ImNo - 1;
        
        if ImNo > 1
            
            Seg_ImAngMid = (Seg_ImAngStartEnd(:,1) + Seg_ImAngStartEnd(:,2))/2;
            
            for I = 1:ImNo-1
                Seg_ImAngDiff(I) = Seg_ImAngMid(I+1) - Seg_ImAngMid(I);
            end
            
        end
        
    end
    
    if strcmp(Plots, 'On')
        
        SegPlot(Image_Diff, N_Sq, N_Cir, N_Seg, R_Cir, xcen, ycen, CirMin, CirMax, Seg_Flux, Seg_Plots, Res_Min2, Center_Type, Im_Cut_Size, PSF_Cut_Size)
       
        FigHandle = figure('Position', [100, 100, 1250, 950], 'Color', [1 1 1]);
        
        plot(Seg_Theta2, abs(Seg_Flux2), 'k-', 'LineWidth',3)
        xlabel('Segment Theta (Degrees)', 'FontSize', 26); ylabel('Absolute Segment Residual Flux', 'FontSize', 26);
        title('Residual Flux profile within Segments', 'FontSize', 26)
        set(gca,'FontSize',26)
        
    end
    
    
    % Final Data on a Candidate
    
    if (Flux_Detect == 1)
        
        %    Cand_Results(1,1:ImNo) = 0;
        
        for I = 1:ImNo
            Cand_Index = Cand_Index + 1;
            Cand_Results(Cand_Index, 1:4) = [Cand,  I, Seg_ImSubtend(I), Seg_ImFlux(I)];
        end
        
       
    end
   
    
end

toc


