clear;
close all

[V, T] = readObj('../mesh_samples/camelhead.obj');

% figure; drawMesh(T, V);

numVertices  = size(V, 1);
numTriangles = size(T, 1);

%% cot laplacian
vecCot = @(x, y) cot(acos(dot(x, y, 2) ./ (vecnorm(x, 2, 2) .* vecnorm(y, 2, 2))));

cot1 = vecCot(V(T(:, 2), :) - V(T(:, 1), :), V(T(:, 3), :) - V(T(:, 1), :));    % cot(∠213)
cot2 = vecCot(V(T(:, 3), :) - V(T(:, 2), :), V(T(:, 1), :) - V(T(:, 2), :));    % cot(∠321)
cot3 = vecCot(V(T(:, 1), :) - V(T(:, 3), :), V(T(:, 2), :) - V(T(:, 3), :));    % cot(∠132)

C = sparse(T(:, [2 3 3 1 1 2]), ...
           T(:, [3 2 1 3 2 1]), ...
           [cot1; cot1; cot2; cot2; cot3; cot3], numVertices, numVertices);
D = sparse(1:numVertices, 1:numVertices, sum(C), numVertices, numVertices);
L = D - C;

dcpL = decomposition(L(2:numVertices, 2:numVertices));

%% tutte
boundary = findBoundary(V, T);
interior = setdiff(1:numVertices, boundary);

u = zeros(numVertices, 1);
u(boundary) = exp(2i* pi * (1:length(boundary)) / length(boundary));

u(interior) = L(interior, interior) \ (-L(interior, boundary) * u(boundary));

% figure; drawMesh(T, u);

%% flip triangles
% areaC = @(z1, z2) real(z1) .* imag(z2) - imag(z1) .* real (z2);
% signedAreas = areaC(u(T(:, 2)) - u(T(:, 1)), u(T(:, 3)) - u(T(:, 1)));
% find(signedAreas < 0)

%% local-global
% local isometric parameterization
x21 = vecnorm(V(T(:, 2), :) - V(T(:, 1), :), 2, 2);
edgeLength31 = vecnorm(V(T(:, 3), :) - V(T(:, 1), :), 2, 2);
cosAng1 = dot(V(T(:, 2), :) - V(T(:, 1), :), V(T(:, 3), :) - V(T(:, 1), :), 2) ./ x21 ./ edgeLength31;
x31 = edgeLength31 .* cosAng1;
x32 = edgeLength31 .* sqrt(1 - cosAng1.^2);

for i = 1:100
    %% local phase
    % Jacobi
    J11 = real( -u(T(:, 1)) + u(T(:, 2)) ) ./ x21;
    J21 = imag( -u(T(:, 1)) + u(T(:, 2)) ) ./ x21;
    J12 = real( (-u(T(:, 1)) + u(T(:, 3))) .* x21 - (-u(T(:, 1)) + u(T(:, 2))) .* x31 ) ./ x21 ./ x32;
    J22 = imag( (-u(T(:, 1)) + u(T(:, 3))) .* x21 - (-u(T(:, 1)) + u(T(:, 2))) .* x31 ) ./ x21 ./ x32;

    % SVD
    [~, ~, theta, phi] = SVD2(J11, J12, J21, J22);
    rotC = exp(1i * (theta + phi));

    %% global phase
    r2c = @(x) x(:, 1) + 1i * x(:, 2);
    zeroF = zeros(numTriangles, 1);
    x1 = zeros(numTriangles, 1);
    x2 = x21;
    x3 = x31 + 1i*x32;
    res = sparse(T, 1, [cot3 .* (x1 - x2) + cot2 .* (x1 - x3), ...
                        cot1 .* (x2 - x3) + cot3 .* (x2 - x1), ...
                        cot2 .* (x3 - x1) + cot1 .* (x3 - x2)] .* rotC, ...
                 numVertices, 1);

    u(2:numVertices) = dcpL \ (res(2:numVertices) - L(2:numVertices, 1) * u(1));
end

figure; drawMesh(T, u);

function [sigma1, sigma2, theta, phi] = SVD2(m11, m12, m21, m22)
    % 2x2 SVD
    % m = rotm(phi) * diag(sigma1, sigma2) * rotm(theta)
    E = (m11 + m22) / 2;
    F = (m11 - m22) / 2;
    G = (m12 + m21) / 2;
    H = (m21 - m12) / 2;
    Q = sqrt(E.^2 + H.^2);
    R = sqrt(F.^2 + G.^2);
    sigma1 = Q + R;
    sigma2 = Q - R;
    a1 = atan2(G, F);
    a2 = atan2(H, E);
    theta = (a2 - a1) / 2;
    phi = (a2 + a1) / 2;
end