function L = laplacianCot(x, f)

nV = size(x, 1);

vecCot = @(x, y) cot(acos(dot(x, y, 2) ./ (vecnorm(x, 2, 2) .* vecnorm(y, 2, 2))));

cot1 = vecCot(x(f(:, 2), :) - x(f(:, 1), :), x(f(:, 3), :) - x(f(:, 1), :));    % cot(∠213)
cot2 = vecCot(x(f(:, 3), :) - x(f(:, 2), :), x(f(:, 1), :) - x(f(:, 2), :));    % cot(∠321)
cot3 = vecCot(x(f(:, 1), :) - x(f(:, 3), :), x(f(:, 2), :) - x(f(:, 3), :));    % cot(∠132)

C = sparse([f(:, 2); f(:, 3); f(:, 3); f(:, 1); f(:, 1); f(:, 2)], ...
           [f(:, 3); f(:, 2); f(:, 1); f(:, 3); f(:, 2); f(:, 1)], ...
           [cot1;    cot1;    cot2;    cot2;    cot3;    cot3], nV, nV);

D = sparse(1:nV, 1:nV, sum(C), nV, nV);

L = (D - C) ./ sum(C)';