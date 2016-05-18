clear all; clc;

% First Define set of circles

N_Cir = 1;

Xc_Cir = 0;
Yc_Cir = 0;

% All Equal Radii distances

for I = 1:N_Cir
R_Cir(I) = 0.5 + (I-1);
end

for I = 1:N_Cir
    Area_Cir(I) = pi*R_Cir(I)^2;
    if (I > 1)
    for J = 2:I
        Area_Cir(I) = Area_Cir(I) - Area_Cir(J-1);
    end
    end
end

%Now Define Square Grid

N_Sq = 3;

% Make N_Sq^2 cell to contain each pixels area fraction

Area_Frac(1:N_Sq, 1:N_Sq, 1:N_Cir) = 0;

funcirc = @(y,r) sqrt(r^2 - y^2)

for I = 1:N_Sq
    for J =1:N_Sq
        
            N_xlim = N_Sq/2. - I;
            N_ylim = N_Sq/2. - J;
        
        for K = 1:N_Cir    
            
           Area = integral(funcirc(, N_xlim, N_xlim+1)
            
        end
    end
end