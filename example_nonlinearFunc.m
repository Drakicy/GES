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

f = @(z) (exp(4 * pi * (z + (sqrt(3) + 1i * sqrt(2)) / 4)) - 1) ./ (exp(4 * pi * (z + (sqrt(3) + 1i * sqrt(2)) / 4)) + 1);

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
        ], sol.PointNum, toc);

    %%

    ax = subplot(2,2,i);

    hold(ax, 'on')
    
    TR = triangulation(sol.DT(:,:), sol.Domain(1,:) + sol.DomainNorm .* sol.DT.Points);
    triplot(TR, '-k');

    for type = [0 -1 1]
        if type == 0
            ind = (sol.CandPoint(:,2) == 0);
        else
            ind = (type * sol.CandPoint(:,2) > 0);
        end

        scatter(...
            real(sol.CandPoint(ind,1)),...
            imag(sol.CandPoint(ind,1)),...
            25, set_color(type+2,:), 'filled');
    end
    
    hold(ax, 'off')

    set(ax, FontSize=fontsize_global);
    xlabel(ax, '$x$', FontSize=fontsize_axes);
    ylabel(ax, '$y$', FontSize=fontsize_axes);
    title(ax, set_title(i), FontSize=fontsize_title);
    
    xlim(x);
    ylim(y);
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
        ], sol.PointNum, toc);

    %%

    ax = subplot(2,2,i);

    hold(ax, 'on')
    
    TR = triangulation(sol.DT(:,:), sol.Domain(1,:) + sol.DomainNorm .* sol.DT.Points);
    triplot(TR, '-k');

    for type = [0 -1 1]
        if type == 0
            ind = (sol.CandPoint(:,2) == 0);
        else
            ind = (type * sol.CandPoint(:,2) > 0);
        end

        scatter(...
            real(sol.CandPoint(ind,1)),...
            imag(sol.CandPoint(ind,1)),...
            25, set_color(type+2,:), 'filled');
    end
    
    hold(ax, 'off')

    set(ax, FontSize=fontsize_global);
    xlabel(ax, '$x$', FontSize=fontsize_axes);
    ylabel(ax, '$y$', FontSize=fontsize_axes);
    title(ax, set_title(i), FontSize=fontsize_title);
    
    xlim(x);
    ylim(y);
end

%%

% exportgraphics(fig, strcat('figures\', mfilename, '2', '.pdf'));
% exportgraphics(fig, strcat('figures\', mfilename, '2','.png'));