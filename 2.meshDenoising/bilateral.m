% function bilateral(V, T , func, numIterations)

numIterations = 2;

%% find the one ring
numTriangles = size(T, 1);
numVertices  = size(V, 1);

V2V = sparse(T, T(:, [2 3 1]), true, numVertices, numVertices);
V2V = V2V & V2V';

V2T = sparse(repmat( (1:numTriangles)', 1, 3 ), T, true, numTriangles, numVertices);
oneRing = V2T(:, T(:, 1)) | V2T(:, T(:, 2)) | V2T(:, T(:, 3));

%% sigma
sigmaS = 1; sigmaC = 1;

%% compute weights for position
c1 = mean(reshape(V(T, 1), numTriangles, 3), 2);
c2 = mean(reshape(V(T, 2), numTriangles, 3), 2);
c3 = mean(reshape(V(T, 3), numTriangles, 3), 2);

cs1 = c1 .* oneRing - c1' .* oneRing;
cs2 = c2 .* oneRing - c2' .* oneRing;
cs3 = c3 .* oneRing - c3' .* oneRing;

wC = spfun( @exp, -(cs1.^2 + cs2.^2 + cs3.^2) / 2 / sigmaC );

%% iteration
for i = 1:numIterations
    fs1 = normals(:, 1) .* oneRing - normals(:, 1)' .* oneRing;
    fs2 = normals(:, 2) .* oneRing - normals(:, 2)' .* oneRing;
    fs3 = normals(:, 3) .* oneRing - normals(:, 3)' .* oneRing;
    wS = spfun( @exp, -(fs1.^2 + fs2.^2 + fs3.^2) / 2 / sigmaS );

    w = wC .* wS;
    w = w ./ sum(w);

    normals = w * normals;
    normals = normals ./ vecnorm(normals, 2, 2);
end