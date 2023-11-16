function result = gaussianCurvature(V, T)
% Gauss-Bonnet Scheme

calcAngles = @(v1, v2) acos(dot(v1, v2, 2) ./ vecnorm(v1, 2, 2) ./ vecnorm(v2, 2, 2));
calcAreas  = @(v1, v2) vecnorm(cross(v1, v2, 2), 2, 2) / 2;

angles = calcAngles( V(T(:, [2 3 1]), :) - V(T, :), V(T(:, [3 1 2]), :) - V(T, :) );
areas  = calcAreas ( V(T(:, [2 3 1]), :) - V(T, :), V(T(:, [3 1 2]), :) - V(T, :) );

result = (2*pi - accumarray(T(:), angles));
% result = (2*pi - accumarray(T(:), angles)) ./ (accumarray(T(:), areas) / 3);