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

f = @(z) (exp(4 * pi * (z + (sqrt(3) + 1i * sqrt(2)) / 4)) - 1) ./ (exp(4 * pi * (z + (sqrt(3) + 1i * sqrt(2)) / 4)) + 1);

x = [-1 1];
y = x;

tol = 1e-3;

%%

point_num_max = [1e2 2e2 4e2 ceil(1 / tol^2)];

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
        sol = GES(f, [x; y]', tol, point_num_max(i));
    end

    %%

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
% exportgraphics(fig, strcat('figures\', mfilename, '1','.png'));

%%

clear('sol')

%%

f = @(z) z + exp(4 * pi * z);

x = [-1 1];
y = x;

tol = 1e-3;

%%

point_num_max = [1e2 2e2 4e2 ceil(1 / tol^2)];

%%

fig = figure(Position=[0 0 1000 1000]);

time_elapsed = 0;

for i = 1:4
    fprintf('Step <strong>%i/4</strong>\n', i);

    tic

    %%

    if exist('sol', 'var')
        sol.PointNumMax = point_num_max(i);
        sol.fitTriang;
    else
        sol = GES(f, [x; y]', tol, point_num_max(i));
    end

    %%

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

% exportgraphics(fig, strcat('figures\', mfilename, '2', '.pdf'));
% exportgraphics(fig, strcat('figures\', mfilename, '2','.png'));