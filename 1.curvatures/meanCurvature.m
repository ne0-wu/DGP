function H = meanCurvature(V, T)
% Paraboloid Fitting

calcNormals = @(v1, v2) cross(v1, v2, 2) ./ vecnorm(cross(v1, v2, 2), 2, 2);

normalsT = calcNormals( V(T(:, [2 3 1]), :) - V(T, :), V(T(:, [3 1 2]), :) - V(T, :) );

normalsV = [accumarray(T(:), normalsT(:, 1)), ...
            accumarray(T(:), normalsT(:, 2)), ...
            accumarray(T(:), normalsT(:, 3))];
normalsV = normalsV ./ vecnorm(normalsV, 2, 2);

V2V = sparse(T, T(:, [2 3 1]), true);
V2V = V2V & V2V';

H = zeros(size(V, 1), 1);

for i = 1:size(V, 1)
    %%
    normal = normalsV(i, :);
    z = [0 0 1];

    rotAxis = cross(normal, z, 2);
    rotAxis = rotAxis ./ vecnorm(rotAxis, 2, 2);

    rotAngle = acos( dot(normal, z, 2) );

    K = [          0 -rotAxis(3)  rotAxis(2);
          rotAxis(3)           0 -rotAxis(1);
         -rotAxis(2)  rotAxis(1)           0];
    R = eye(3) + sin(rotAngle) * K + (1-cos(rotAngle)) * K^2;

    %%
    nbV = (V(V2V(:, i), :) - V(i, :)) * R';
    coe = [nbV(:, 1) .^ 2, nbV(:, 1) .* nbV(:, 2), nbV(:, 2) .^ 2] \ nbV(:, 3);
    H(i) = coe(1) + coe(3);
end

% visualizeCurvature(V, T, H);