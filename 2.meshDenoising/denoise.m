clear
close all

[V, T] = readObj('../mesh_samples/cathead.obj');

%% original mesh's normals
origNormals = cross(V(T(:, 2), :) - V(T(:, 1), :), V(T(:, 3), :) - V(T(:, 1), :));
origNormals = origNormals ./ vecnorm(origNormals, 2, 2);

%% add noise to the mesh

%% compute normals
normals = cross(V(T(:, 2), :) - V(T(:, 1), :), V(T(:, 3), :) - V(T(:, 1), :));
normals = normals ./ vecnorm(normals, 2, 2);

%% apply bilateral filter
bilateral

%% fit the normals

% T2 = T(:, [2 3 1]);
% A = sparse([repmat(1:3*numTriangles, 1, 6)], ...
%            [T, T2, T+numVertices, T2+numVertices, T+2*numVertices, T2+2*numVertices], ...
%            [normals(:, [1 1 1]), -normals(:, [1 1 1]), normals(:, [2 2 2]), -normals(:, [2 2 2]), normals(:, [3 3 3]), -normals(:, [3 3 3])]);
% 
% boundary = findBoundary(V, T);
% lambda = 1e3;
% B1 = lambda * sparse(1:length(boundary), boundary, ones(1, length(boundary)), length(boundary), numVertices*3);
% B2 = lambda * sparse(1:length(boundary), numVertices+boundary, ones(1, length(boundary)), length(boundary), numVertices*3);
% B3 = lambda * sparse(1:length(boundary), numVertices*2+boundary, ones(1, length(boundary)), length(boundary), numVertices*3);
% 
% V2 = [A; B1; B2; B3] \ [zeros(3*numTriangles, 1); lambda * reshape(V(boundary, :), length(boundary)*3, 1)];
% V2 = reshape(V2, numVertices, 3);
% 
% trimesh(T, V2(:, 1), V2(:, 2), V2(:, 3))

c = [c1 c2 c3];

for count = 1:5
    for i = 1:numVertices
        nbTriangles = find(V2T(:, i));
        numNbTriangles = length(nbTriangles);
        for j = nbTriangles'
            V(i, :) = V(i, :) + (c(j, :) - V(i, :)) * normals(j, :)' * normals(j, :) / numNbTriangles;
        end
    end
end

figure
trimesh(T, V(:, 2), V(:, 1), V(:, 3))
axis equal