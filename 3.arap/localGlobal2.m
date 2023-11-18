clear
close all

[V, T] = readObj('../mesh_samples/camelhead.obj');

% figure; drawMesh(V, F);

numVertices = size(V, 1);
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

%% tutte
boundary = findBoundary(V, T);
interior = setdiff(1:numVertices, boundary);

u = zeros(numVertices, 1);
u(boundary) = exp(2i* pi * (1:length(boundary)) / length(boundary));

u(interior) = L(interior, interior) \ (-L(interior, boundary) * u(boundary));

%% 平面化，求角度
% xs = zeros(numTriangles, 3);  % 网格中三角形放到平面上
% thetas = zeros(numTriangles, 3);  % 三角形中角的弧度
% 
% cross2 = @(a, b) a(1) * b(2) - a(2) * b(1);    % 二维叉积，定向用
% 
% for f = 1:numTriangles
%     vertex = V(T(f, :), :);
%     edge = vertex([3 1 2], :) - vertex([2 3 1], :);    % edges(i, :) 为顶点 i 所对的边
%     thetas(f, :) = acos(dot(edge([3 1 2], :), -edge([2 3 1], :), 2) ./ (vecnorm(edge([3 1 2], :), 2, 2) .* vecnorm(edge([2 3 1], :), 2, 2))); % thetas(t, i) 为顶点 i 所对应的角
%     xs(f, 1) = vecnorm(edge(3, :));
%     if cross2(u(T(f, 2), :) - u(T(f, 1), :), u(T(f, 3), :) - u(T(f, 1), :)) > 0 % 根据 Tutte 参数化中的三角形方向决定铺到平面上的三角形方向
%         xs(f, 2:3) = vecnorm(edge(2, :)) * [cos(thetas(f, 1)), sin(thetas(f, 1))];
%     else
%         xs(f, 2:3) = vecnorm(edge(2, :)) * [cos(thetas(f, 1)), sin(-thetas(f, 1))]; % 事实上对于定向好的网格这是不必要的
%     end
% end
% cots = cot(thetas);

% figure; drawMesh(T, u);

x21 = vecnorm(V(T(:, 2), :) - V(T(:, 1), :), 2, 2);
edgeLength31 = vecnorm(V(T(:, 3), :) - V(T(:, 1), :), 2, 2);
cosAng1 = dot(V(T(:, 2), :) - V(T(:, 1), :), V(T(:, 3), :) - V(T(:, 1), :), 2) ./ x21 ./ edgeLength31;
x31 = edgeLength31 .* cosAng1;
x32 = edgeLength31 .* sqrt(1 - cosAng1.^2);

xs = [x21 x31 x32];


decomp = decomposition(L(2:numVertices, 2:numVertices));

%% local-global
myr2c = @(x) x(:, 1) + 1i * x(:, 2);    % Real vectors to complex numbers
bool2pn = @(b) (b - 0.5) * 2;   % 1 => 1, 0 -> -1

N = 10;

u = [real(u) imag(u)]

for it = 1:N
    J11 = (-u(T(:, 1), 1) + u(T(:, 2), 1)) ./ xs(:, 1);
    J12 = (-u(T(:, 1), 2) + u(T(:, 2), 2)) ./ xs(:, 1);
    J21 = (u(T(:, 3), 1) .* xs(:, 1) - u(T(:, 2), 1) .* xs(:, 2) + u(T(:, 1), 1) .* (-xs(:, 1) + xs(:, 2))) ./ (xs(:, 1) .* xs(:, 3));
    J22 = (u(T(:, 3), 2) .* xs(:, 1) - u(T(:, 2), 2) .* xs(:, 2) + u(T(:, 1), 2) .* (-xs(:, 1) + xs(:, 2))) ./ (xs(:, 1) .* xs(:, 3));
    
    %% faster SVD from https://lucidar.me/en/mathematics/singular-value-decomposition-of-a-2x2-matrix/
    as = atan2(2 * J11 .* J21 + 2 * J12 .* J22, J11.^2 + J12.^2 - J21.^2 - J22.^2) / 2;
    bs = atan2(2 * J11 .* J12 + 2 * J21 .* J22, J11.^2 - J12.^2 + J21.^2 - J22.^2) / 2;

    signcosb = ((cos(as) .* J11 + sin(as) .* J21) .* cos(bs) + (cos(as) .* J12 + sin(as) .* J22) .* sin(bs)) > 0;
    signsinb = ((sin(as) .* J11 - cos(as) .* J21) .* sin(bs) + (-sin(as) .* J12 + cos(as) .* J22) .* cos(bs)) > 0;
    signdetJ = J11 .* J22 - J12 .* J21 > 0;
    bs = atan2(bool2pn(signdetJ) .* bool2pn(signsinb) .* sin(bs), bool2pn(signcosb) .* cos(bs));

    %%
    Ls = exp(1i * (bs - as));
    
    % 这样创建 rhs 是为了利用 sparse 会把相同指标的项加起来的特性
    rhs = full(sparse(T, 1, [cot3 .* myr2c(-[xs(:, 1), zeros(numTriangles, 1)]) .* Ls + cot2 .* myr2c(-xs(:, 2:3)) .* Ls; ...
                           cot1 .* myr2c([xs(:, 1), zeros(numTriangles, 1)] - xs(:, 2:3)) .* Ls + cot3 .* myr2c([xs(:, 1), zeros(numTriangles, 1)]) .* Ls; ...
                           cot2 .* myr2c(xs(:, 2:3)) .* Ls + cot1 .* myr2c(xs(:, 2:3) - [xs(:, 1), zeros(numTriangles, 1)]) .* Ls], ...
                       numVertices, 1));

%     u = [laplacianC; 1, zeros(1, nV-1)] \ [real(rhs), imag(rhs); u(1, :)];
%     u = r \ (r' \ [real(rhs), imag(rhs)]);
%     u(2:nV, :) = laplacianC(2:nV, 2:nV) \ ([real(rhs(2:nV)), imag(rhs(2:nV))] - laplacianC(2:nV, 1) * u(1, :));
    u(2:numVertices, :) = decomp \ ([real(rhs(2:numVertices)), imag(rhs(2:numVertices))] - L(2:numVertices, 1) * u(1, :));
end

figure; drawMesh(T, u);