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

fn = @(z,n) double(psi(n,sym(z)));

x = [-3.75 3.75];
y = [-1 1];

tol = 1e-3;
point_num_max = 0;

%%

figure(Position=[0 0 1000 500]);

for i = 1:2
    fprintf('Step <strong>%i/2</strong>\n', i);

    tic

    %%

    sol = GES(@(z) fn(z, i - 1), [x; y]', tol, PointNumMax=point_num_max, PropMax=1);

    %%

    fprintf(...
        [
            'Points number: <strong>%i</strong>\n'...
            'Overall time: <strong>%f sec</strong>\n'
        ], sol.PointNum, toc);

    %%

    ax = subplot(1,2,i);

    sol.visTriang([0 -1 1]);
    
    title(ax, set_title(i), FontSize=22);
end

%%

% exportgraphics(fig, strcat('figures\', mfilename, '.pdf'));
% exportgraphics(fig, strcat('figures\', mfilename, '.png'));