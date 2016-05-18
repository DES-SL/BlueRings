 clear all; clc;

% First Define set of circles

N_Cir = 10;

Xc_Cir = 0;
Yc_Cir = 0;

% All Equal Radii distances

R_Cir(1) = 0;

for I = 2:N_Cir+1
R_Cir(I) = 0.5 + (I-1);
end

for I = 1:N_Cir
    Area_Cir(I) = pi*R_Cir(I+1)^2;
    if (I > 1)
    for J = 2:I
        Area_Cir(I) = Area_Cir(I) - Area_Cir(J-1)
    end
    end
end

%Now Define Square Grid

N_Sq = 21;

% Make N_Sq^2 cell to contain each pixels area fraction

Area_Frac(1:N_Sq, 1:N_Sq, 1:N_Cir) = 0;

N_Tot = 5000;

for I = 1:N_Sq
    for J =1:N_Sq
        for K = 2:N_Cir+1        
            
            N_Count = 0;
            N_Shiftx = N_Sq/2. - I + 1;
            N_Shifty = N_Sq/2. - J + 1;
            
            %Now do Monte Carlo
            for N = 1:N_Tot
                randx = rand - N_Shiftx;
                randy = rand - N_Shifty;
                
                % If inside circle count 1
                if (randx^2 + randy^2) < R_Cir(K)^2
                    if (randx^2 + randy^2) > R_Cir(K-1)^2
                    N_Count = N_Count+1;
                    end
                end
            end
              
            Area_Frac(I,J,K-1) = N_Count/N_Tot

        end
    end
end