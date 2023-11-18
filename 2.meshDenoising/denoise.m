clear
close all

[V, T] = readObj('../mesh_samples/ball.obj');

numBilateralIterations = 2000;
numFitIterations = 20;
sigmaS = 1;

%% original mesh's normals
origNormals = cross(V(T(:, 2), :) - V(T(:, 1), :), V(T(:, 3), :) - V(T(:, 1), :));
origNormals = origNormals ./ vecnorm(origNormals, 2, 2);

%% add noise to the mesh

V = V + randn(size(V)) * mean(vecnorm(V(T, :) - V(T(:, [2 3 1]), :), 2, 2)) / 5;

figure
trimesh(T, V(:, 1), V(:, 2), V(:, 3));
axis equal

%% compute normals
normals = cross(V(T(:, 2), :) - V(T(:, 1), :), V(T(:, 3), :) - V(T(:, 1), :));
normals = normals ./ vecnorm(normals, 2, 2);

%% apply bilateral filter to normals
bilateral

%% fit the normals

c = [c1 c2 c3];

for count = 1:numFitIterations
    for i = 1:numVertices
        nbTriangles = find(V2T(:, i));
        numNbTriangles = length(nbTriangles);
        for j = nbTriangles'
            V(i, :) = V(i, :) + (c(j, :) - V(i, :)) * normals(j, :)' * normals(j, :) / numNbTriangles;
        end
    end
end

figure
trimesh(T, V(:, 1), V(:, 2), V(:, 3))
axis equal