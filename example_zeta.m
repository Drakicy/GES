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

f = @(z) zeta(z);

x = [-9.75 9.75];
y = 2 * x;

tol = 1e-3;

%%

point_num_max = [1e2 2e2 4e2 0];

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
        sol = GES(f, [x; y]', tol, PointNumMax=point_num_max(i));
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

% exportgraphics(fig, strcat('figures\', mfilename, '.pdf'));
% exportgraphics(fig, strcat('figures\', mfilename, '.png'));