function fo = nan_interpolate(go)

% Interpolate NaN values
%
% Input
% go: data with NaN data (1xN)
%
% Output
% fo: data without NaN data (1xN)

x = find(~isnan(go));
y = go(x);
xi = 1 : max(length(go));
fo = interp1(x, y, xi);
fo = transpose(fo);

end