function [Area_Frac, Area_Cir, Cir_Cross, R_Cir, R_Pix_Max, R_Pix_Min, Seg_Frac, Seg_Grid, Seg_Theta] = MC_Areas(N_Cir, R_Cir_Space, Xc_Cir, Yc_Cir, N_Sq, N_Seg, N_Tot)

% First calc the radii of all circles

for I = 1:N_Cir
    R_Cir(I) = R_Cir_Space*I;
end

% You compair different annuli for area fractions, this moves the array so
% that circle 1 compares to a R = 0

R_Cir_Comp(1) = 0;
R_Cir_Comp(2:N_Cir+1) = R_Cir(:);

%Now calc their annuli areas, but calculating each circles total area and
%subtracting off all other circle areas

for I = 1:N_Cir
    Area_Cir(I) = pi*R_Cir(I)^2;
    if (I > 1)
        for J = 2:I
            Area_Cir(I) = Area_Cir(I) - Area_Cir(J-1);
        end
    end
end

% Make N_Sq^2 cell to contain each pixels area fraction

Area_Frac(1:N_Sq, 1:N_Sq, 1:N_Cir, 1:2) = 0;
Seg_Frac(1:N_Sq, 1:N_Sq, 1:2) = 0;
R_Pix_Grid(1:N_Sq, 1:N_Sq, 1:2) = 0;

% Center Type 1, where center is at the middle of a pixel

for I = 1:N_Sq
    for J =1:N_Sq
        for K = 1:N_Cir
            
            N_Count = 0;
            N_Shiftx = N_Sq/2. - I + 1;
            N_Shifty = N_Sq/2. - J + 1;
            
            %Now do Monte Carlo
            for N = 1:N_Tot
                
                randx = rand - N_Shiftx;
                randy = rand - N_Shifty;
                
                % If inside circle count 1
                if ((randx+Xc_Cir)^2 + (randy+Yc_Cir)^2) < R_Cir_Comp(K+1)^2
                    if ((randx+Xc_Cir)^2 + (randy+Yc_Cir)^2) > R_Cir_Comp(K)^2
                        N_Count = N_Count+1;
                    end
                end
            end
            
            Area_Frac(I,J,K,1) = N_Count/N_Tot;
            
            if Area_Frac(I,J,K,1) > 0 && Area_Frac(I,J,K,1) < 1
                
                Cir_Cross(I,J,1) = K;
                
            end
            
        end
    end
end

for I = 1:N_Sq
    for J = 1:N_Sq
        
        N_Count = 0;
        N_Shiftx = N_Sq/2. - I + 1;
        N_Shifty = N_Sq/2. - J + 1;
        
        %Calc the max,min distance each pixels corners are from the center
        
        r1 = sqrt((N_Shiftx - 0.5 + 0.5)^2 + (N_Shifty - 0.5 + 0.5)^2);
        r2 = sqrt((N_Shiftx - 0.5 - 0.5)^2 + (N_Shifty - 0.5 + 0.5)^2);
        r3 = sqrt((N_Shiftx - 0.5 + 0.5)^2 + (N_Shifty - 0.5 - 0.5)^2);
        r4 = sqrt((N_Shiftx - 0.5 - 0.5)^2 + (N_Shifty - 0.5 - 0.5)^2);
        
        R_Pix_Max(I,J,1) = max([r1, r2, r3, r4]);
        R_Pix_Min(I,J,1) = min([r1, r2, r3, r4]);
        
        
        % If it is within, calculate its segment
        theta = mod(atan2((0.5 - N_Shifty),(0.5 - N_Shiftx)),2*pi);
        Seg_No = ceil((theta / (2*pi))*N_Seg);
        
        % Thetas of two neighbouring segments
        Seg_theta(1) = (Seg_No-1)*((2*pi)/N_Seg);
        Seg_theta(2) = (Seg_No)*((2*pi)/N_Seg);
        
        %Compare areas with two neighbour segment lines
        
        for N = 1:N_Tot
            
            randx = rand - N_Shiftx;
            randy = rand - N_Shifty;
            
            rand_theta = mod(atan2(randy,randx),2*pi);
            
            
            % If inside Arc count 1
            if (rand_theta > Seg_theta(1))
                if (rand_theta < Seg_theta(2))
                    N_Count = N_Count+1;
                end
            end
        end
        
        Seg_Frac(I,J,1) = N_Count/N_Tot;
        
        % Seg Grid says which segment each pixel belongs to
        
        Seg_Grid(I,J,1) = Seg_No;
        
        if (Seg_Frac(I,J,1) < 1)
            
            if (((theta / (2*pi))*N_Seg) - Seg_No) < -0.5
                
                Seg_Grid(I,J,1) = Seg_No-0.5;
                
            elseif (((theta / (2*pi))*N_Seg) - Seg_No) > -0.5
                Seg_Grid(I,J,1) = Seg_No+0.5;
            end
            
        end
        
    end
end

cenrow = (N_Sq/2)+0.5;

Seg_Frac(cenrow, cenrow) = 1/N_Seg;
Seg_Frac(cenrow, 1:cenrow) = fliplr(Seg_Frac(cenrow, cenrow:N_Sq));
Seg_Frac(cenrow:N_Sq, cenrow) = flipud(Seg_Frac(1:cenrow, cenrow));

% Center Type 2, is a grid where the 'centre' is in the middle of 4 pixels,
% thus the grid must be shifted appropriately

for I = 1:N_Sq
    for J =1:N_Sq
        for K = 1:N_Cir
            
            N_Count = 0;
            N_Shiftx = N_Sq/2. - I + 0.5;
            N_Shifty = N_Sq/2. - J + 0.5;
            
            %Now do Monte Carlo
            for N = 1:N_Tot
                randx = rand - N_Shiftx;
                randy = rand - N_Shifty;
                
                % If inside circle count 1
                if ((randx+Xc_Cir)^2 + (randy+Yc_Cir)^2) < R_Cir_Comp(K+1)^2
                    if ((randx+Xc_Cir)^2 + (randy+Yc_Cir)^2) > R_Cir_Comp(K)^2
                        N_Count = N_Count+1;
                    end
                end
            end
            
            Area_Frac(I,J,K,2) = N_Count/N_Tot;
            
            if Area_Frac(I,J,K,2) > 0 && Area_Frac(I,J,K,2) < 1
                
                Cir_Cross(I,J,2) = K;
                
            end
            
        end
    end
end

for I = 1:N_Sq
    for J = 1:N_Sq
        
        N_Count = 0;
        N_Shiftx = N_Sq/2. - I + 0.5;
        N_Shifty = N_Sq/2. - J + 0.5;
        
        %Calc the max,min distance each pixels corners are from the center
        
        r1 = sqrt((N_Shiftx - 0.5 + 0.5)^2 + (N_Shifty - 0.5 + 0.5)^2);
        r2 = sqrt((N_Shiftx - 0.5 - 0.5)^2 + (N_Shifty - 0.5 + 0.5)^2);
        r3 = sqrt((N_Shiftx - 0.5 + 0.5)^2 + (N_Shifty - 0.5 - 0.5)^2);
        r4 = sqrt((N_Shiftx - 0.5 - 0.5)^2 + (N_Shifty - 0.5 - 0.5)^2);
        
        R_Pix_Max(I,J,2) = max([r1, r2, r3, r4]);
        R_Pix_Min(I,J,2) = min([r1, r2, r3, r4]);
        
        
        % Is square pixel within the circular annuli?
        %    R_Pix_Grid(I,J,2) = sqrt((N_Shiftx - 0.5)^2 + (N_Shifty - 0.5)^2);
        
        % If it is within, calculate its segment
        theta = mod(atan2((0.5 - N_Shifty),(0.5 - N_Shiftx)),2*pi);
        Seg_No = ceil((theta / (2*pi))*N_Seg);
        
        % Thetas of two neighbouring segments
        Seg_theta(1) = (Seg_No-1)*((2*pi)/N_Seg);
        Seg_theta(2) = (Seg_No)*((2*pi)/N_Seg);
        
        %Compare areas with two neighbour segment lines
        
        for N = 1:N_Tot
            
            randx = rand - N_Shiftx;
            randy = rand - N_Shifty;
            
            rand_theta = mod(atan2(randy,randx),2*pi);
            
            
            % If inside Arc count 1
            if (rand_theta > Seg_theta(1))
                if (rand_theta < Seg_theta(2))
                    N_Count = N_Count+1;
                end
            end
        end
        
        Seg_Frac(I,J,2) = N_Count/N_Tot;
        
        % Seg Grid says which segment each pixel belongs to
        
        Seg_Grid(I,J,2) = Seg_No;
        
        if (Seg_Frac(I,J,2) < 1)
            
            if (((theta / (2*pi))*N_Seg) - Seg_No) < -0.5
                Seg_Grid(I,J,2) = Seg_No-0.5;
            elseif (((theta / (2*pi))*N_Seg) - Seg_No) > -0.5
                Seg_Grid(I,J,2) = Seg_No+0.5;
            end
            
        end
        
    end
end

for I = 1:N_Seg
    Seg_Theta(I,1) = (I-1)*((2*pi)/N_Seg) + ((2*pi)/N_Seg)/2;
end

for I = 1:N_Sq
    for J = 1:N_Sq
        for K = 1:2
            
            if (Seg_Grid(I,J,K) > N_Seg)
                Seg_Grid(I,J,K) = N_Seg
            end
            
        end
    end
end

%Area_Frac = Area_Frac(1:N_Sq-1, 1:N_Sq-1, :);

end