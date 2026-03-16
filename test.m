
target_idx = 5;

f = 1575.42e6;
C = 299792458;

a = diff(pr1(:, target_idx));
b = -dop1(1:end-1, target_idx) * C / f;
c = diff(ph1(:, target_idx)) * (C/f);

idx = 1:length(a);

figure(3);
clf;

hold on;
scatter(idx, a, 'Marker','o', 'MarkerEdgeColor','blue');
scatter(idx, c, 'Marker', '+', 'MarkerEdgeColor','red');
scatter(idx, b);
