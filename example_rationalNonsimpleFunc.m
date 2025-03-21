clear
set(0, ...
        defaultTextInterpreter='latex', ...
        defaultAxesTickLabelInterpreter='latex', ...
        defaultAxesFontSize=18, ...
        defaultAxesLabelFontSize=1.4, ...
        defaultAxesTitleFontSize=1.4, ...
        defaultLegendInterpreter='latex', ...
        defaultLegendFontSize=16, ...
        defaultLineLineWidth=1, ...
        defaultLineMarkerSize=12 ...
    );

color = ...
    [
    0.8500 0.3250 0.0980;
    0 0 0;
    0 0.4470 0.7410
    ];

%%

f = @(z,a) (z - a).^3 ./ (z + a).^2;

x = [-1 1];
y = x;

a = (sqrt(3) + 1i * sqrt(2)) / 4;

%%

point_num_max = [1e1 5e1 1e2 0];

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
        sol = GES( ...
                @(z) f(z,a), ...
                [x; y]', ...
                PointNumMax=point_num_max(i) ...
            );
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

    sol.visTriang([0 -1 1]);
    
    title("(" + char('a' + (i - 1)) + ")");
end

%%

% exportgraphics(fig, strcat('figures\', mfilename, '.pdf'));
% exportgraphics(fig, strcat('figures\', mfilename, '.png'));