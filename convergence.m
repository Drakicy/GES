clear
set(0, defaultAxesTickLabelInterpreter='latex');
set(0, defaulttextinterpreter='latex');
set(0, defaultLegendInterpreter='latex');

fontsize_global = 18;
fontsize_axes = 25;
fontsize_title = 22;

set_title = ["(a)" "(b)" "(c)" "(d)"];
set_color = ...
    [
    0.8500 0.3250 0.0980;
    0 0 0;
    0 0.4470 0.7410
    ];

%%

f = @(z,a) z - a;

x = [-1 1];
y = x;

a = (sqrt(3) + 1i * sqrt(2)) / 4;
tol = 1e-14;

%%

n = 100;
point_num_max = ceil(linspace(1e2, 5e2, n));

approx_error = zeros(n,1);

for i = 1:n
    fprintf('Step <strong>%i/%i</strong>\n', i, n);

    if exist('sol', 'var')
        sol.PointNumMax = point_num_max(i);
        sol.fitTriang;
    else
        sol = GES(@(z) f(z,a), [x; y]', tol, point_num_max(i));
    end

    approx_error(i) = abs(mean(sol.CandPoint(:,1)) - a);
end

%%

fig = figure(Position=[0 0 500 500]);

error_fit = fit(point_num_max', log(approx_error), 'poly1');
coef = coeffvalues(error_fit);

hold on

data = plot(...
        point_num_max,...
        approx_error,...
        '.b'...
    );
data_fit = plot(...
        point_num_max,...
        exp(error_fit(point_num_max)),...
        '-r',...
        LineWidth=2 ...
    );

hold off

set(gca,...
    FontSize=fontsize_global,...
    YScale='log'...
    );
xlabel(gca, '$N_{tr}$', FontSize=fontsize_axes);
ylabel(gca, '$|z_i - a|$', FontSize=fontsize_axes);
legend(...
        [data data_fit],...
        'Approximation error',...
        sprintf('$%.2f \\exp(%.2f N_{tr})$', exp(coef(2)), coef(1))...
    )