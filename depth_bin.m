function [u_bin, dep_bin] = depth_bin(u, dep, max_dep, dz)

% Binning data in water depth
%
% Input:
% u:        Data for binning
% dep:      Water depth from surface
% max_dep:  Maximum depth
% dz:       Depth interval for binning
%
% Output:
% u_bin:     Binned velocity component u
% v_bin:     Binned velocity component v
% u_dep_ave: Depth-averaged velocity component u
% v_dep_ave: Depth-averaged velocity component v
% w_fac:     Weight factor to calculate depth-averaged velocity

%           ------------------ h_1  | z direction
% 1-st bin    h_1 <= z < h_2        v
%           ------------------ h_2            
% 2-nd bin    h_2 <= z < h_3
%           ------------------ h_3

% Number of bins

n_bin = round(max_dep / dz);

% Water depth of each bin

dep_bin = zeros(n_bin, 1);
for z = 1 : size(dep_bin, 1)
    dep_bin(z) = 0.5*dz + (z-1)*dz;
end

% --- Binning data

u_bin = 0;
for bin = 1 : n_bin  % Loop for bin
    
    h_i = dz * (bin - 1) - 1.e-06;
    
    count = 0;
    u_bar = 0;
    for i = 1 : size(u, 1)
        if dep(i) >= h_i && dep(i) < h_i+dz
            count = count + 1;
            u_bar = u_bar + u(i);
        end
    end
    if count > 0
        u_bar = u_bar / count;
    else
        u_bar = NaN;
    end
    u_bin(bin,1) = u_bar;
end

% --- Interpolate NaN values

go = u_bin;
interp = nan_interpolate(go);
u_bin = interp;

end