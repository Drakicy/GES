clear
set(0, defaultAxesTickLabelInterpreter='latex');
set(0, defaulttextinterpreter='latex');
set(0, defaultLegendInterpreter='latex');
set_title = ["(a)" "(b)" "(c)" "(d)"];
set_color = ...
    [
    0.8500 0.3250 0.0980;
    0 0 0;
    0 0.4470 0.7410
    ];

%%

f = @(z,a) (z - a) ./ (z + a);

x = [-1 1];
y = x;

a = (sqrt(3) + 1i * sqrt(2)) / 4;
tol = 1e-3;

%%

point_num_max = [5e1 1e2 2e2 0];

%%

figure(Position=[0 0 1000 1000]);

time_elapsed = 0;

for i = 1:4
    fprintf('Step <strong>%i/4</strong>\n', i);

    tic

    %%

    if exist('sol', 'var')
        sol.PointNumMax = point_num_max(i);
        sol.fitTriang;
    else
        sol = GES(@(z) f(z,a), [x; y]', tol, PointNumMax=point_num_max(i));
    end

    time_elapsed = time_elapsed + toc;

    fprintf(...
        [
            'Points number: <strong>%i</strong>\n'...
            'Overall time: <strong>%f sec</strong>\n'
        ], sol.PointNum, time_elapsed);

    %%

    ax = subplot(2,2,i);

    sol.visTriang;
    
    title(ax, set_title(i), FontSize=22);
end

%%

% exportgraphics(fig, strcat('figures\', mfilename, '1', '.pdf'));
% exportgraphics(fig, strcat('figures\', mfilename, '1', '.png'));

%%

a = (sqrt(3) + 1i * sqrt(2)) / 4 * [1e-1 1e-2];
point_num_max = 0;

%%

fig = figure(Position=[0 0 1000 500]);

for i = 1:2
    fprintf('Step <strong>%i/2</strong>\n', i);

    tic

    %%

    sol = GES(@(z) f(z,a(i)), [x; y]', tol, PointNumMax=point_num_max);

    fprintf(...
        [
            'Points number: <strong>%i</strong>\n'...
            'Overall time: <strong>%f sec</strong>\n'
        ], sol.PointNum, toc);

    %%

    ax = subplot(1,2,i);

    sol.visTriang;
    
    title(ax, set_title(i), FontSize=22);
end

%%

% exportgraphics(fig, strcat('figures\', mfilename, '2', '.pdf'));
% exportgraphics(fig, strcat('figures\', mfilename, '2', '.png'));