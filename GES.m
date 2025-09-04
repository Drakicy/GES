classdef GES < handle
    properties (SetAccess=private)
        %% Domain: domain boundaries
        Domain (2,2)
        %   matrix

        %% DomainNorm: domain normalization
        DomainNorm (1,1)
        %   scalar

        %% Func: function
        Func function_handle
        %   function handle

        %% DT: Delaunay triangulation
        DT = delaunayTriangulation()
        %   Delaunay triangulation class

        %% FuncEval: function evaluations
        FuncEval (:,1)
        %   vector

        %% PointNum: number of triangulation points
        PointNum (1,1)
        %   scalar

        %% CandPoint: candidate points
        CandPoint (:,:)
        %   matrix

        %% CandRegion: candidate regions
        CandRegion (:,1)
        %   cell array
    end

    properties
        %% Tol: approximation relative tolerance level
        Tol (1,1) {mustBePositive, mustBeLessThanOrEqual(Tol, 1)} = 1
        %   scalar

        %% PointNumMax: maximum number of triangulation points
        PointNumMax (1,1) {mustBeInteger, mustBeNonnegative} = 0
        %   scalar

        %% PointNumMin: minimum number of triangulation points
        PointNumMin (1,1) {mustBeInteger, mustBeNonnegative} = 9
        %   scalar

        %% FlowMax: maximum of gradient flow
        FlowMax (1,1) {mustBePositive} = Inf
        %   scalar

        %% BatchSize: batch size
        BatchSize (1,1) {mustBeInteger, mustBeNonnegative} = 0
        %   scalar

        %% Display: output flag
        Display {mustBeMember(Display, ["off" "on"])} = 'off'
        %   either 'off' or 'on'
    end

    properties (Access=private, Hidden)
        %% SimplexGrad: simplex gradient
        SimplexGrad (:,6) = zeros(0,6)
        %   matrix

        %% PointStorage: leftover points storage
        PointStorage (:,:)
        %   matrix

        %% StartTime: start time
        StartTime (1,1)
        %   scalar
    end

    methods
        function obj = GES(func, dom, opt)
            %% Global Equation Solver (GES)
            %
            % INPUT
            %
            %   func: handle of a function f = f(z) = f(x + iy)
            %
            %   dom: domain boundaries,
            %        2x2 matrix with ascending elements in columns
            %
            %   (optional) Tol: approximation relative tolerance level,
            %                   positive scalar with value less than 1
            %   default: 1 / 2^10
            %
            %   (optional) PointNumMax: maximum number of triangulation points,
            %                           positive integer
            %   default: 0 (no limit)
            %
            %   (optional) PointNumMin: minimum number of triangulation points,
            %                           positive integer
            %   default: empty (2 * (1 + relation between domain sides) + 5)
            %
            %   (optional) FlowMax: maximum of gradient flow,
            %                       positive scalar
            %   default: Inf
            %
            %   (optional) AddPoint: additional initial points,
            %                        complex vector
            %   default: empty
            %
            %   (optional) BatchSize: batch size,
            %                         nonnegative integer
            %   default: 0 (no batching)
            %
            %   (optional) Display: output flag,
            %                       either 'off' or 'on'
            %   default: 'off'
            %
            % OUTPUT
            %
            %   obj: GES class
        
            arguments
                func function_handle
                dom (2,2) {mustBeReal}
                opt.Tol (1,1) {mustBePositive, mustBeLessThan(opt.Tol,1)} = 1 / 2^10
                opt.PointNumMax (1,1) {mustBeInteger, mustBeNonnegative} = 0
                opt.PointNumMin {mustBeScalarOrEmpty, mustBeInteger, mustBeNonnegative} = []
                opt.FlowMax (1,1) {mustBePositive} = Inf
                opt.AddPoint (:,1) = []
                opt.BatchSize (1,1) {mustBeInteger, mustBeNonnegative} = 0
                opt.Display {mustBeMember(opt.Display, ["off" "on"])} = 'off'
            end

            %% Set class properties and perform argument validation

            % Check whether domain boundaries are sorted
            assert( ...
                    issorted(dom, 1), ...
                    'Data:notSorted', ...
                    "Columns of 'dom' must be sorted." ...
                );

            % Check whether domain boundaries are different
            assert( ...
                    all(diff(dom,1) ~= 0), ...
                    'Data:notUnique', ...
                    "Row elements of 'dom' must be different." ...
                );

            % Set domain boundaries
            obj.Domain = dom;

            % Set domain normalization
            side = diff(dom, 1, 1);
            obj.DomainNorm = min(side);
            side_rel = side / obj.DomainNorm;

            % Create initial mesh
            [initial_point_x, initial_point_y] = ...
                ndgrid( ...
                        [0 side_rel(1)], ...
                        [0 side_rel(2)] ...
                    );

            % Check whether additional initial points are inside the domain
            assert( ...
                    all((real(opt.AddPoint) >= dom(1,1)) & (real(opt.AddPoint) <= dom(2,1)) & (imag(opt.AddPoint) >= dom(1,2)) & (imag(opt.AddPoint) <= dom(2,2))), ...
                    'Data:outOfDomain', ...
                    "Elements of 'AddPoint' must be inside the domain." ...
                );

            % Set initial triangulation points
            obj.PointStorage = ...
                unique( ...
                    [
                        initial_point_x(:) initial_point_y(:)
                        ([real(opt.AddPoint) imag(opt.AddPoint)] - obj.Domain(1,:)) / obj.DomainNorm
                    ] ...
                , 'rows');

            % Set approximation relative tolerance level
            obj.Tol = opt.Tol;

            % Set function
            obj.Func = @(z) func((obj.Domain(1,:) + obj.DomainNorm .* z) * [1; 1i]);

            % Set maximum number of triangulation points
            obj.PointNumMax = opt.PointNumMax;

            % Set minimum number of triangulation points
            if isempty(opt.PointNumMin)
                obj.PointNumMin = 2 * (1 + floor(prod(side_rel))) + 5;
            else
                obj.PointNumMin = opt.PointNumMin;
            end

            % Set maximum of gradient flow
            obj.FlowMax = opt.FlowMax;

            % Set batch size
            obj.BatchSize = opt.BatchSize;

            % Set output flag
            obj.Display = opt.Display;

            %% Perform triangulation fitting

            % Initiate triangulation fitting
            obj.fitTriang;
        end

        function fitTriang(obj)
            %% fitTriang: triangulation fitting
            %
            % INPUT
            %
            % OUTPUT
            %

            % Set start time
            obj.StartTime = tic;

            %% Complete previous run

            if ~isempty(obj.PointStorage)
                % Add leftover triangulation points
                obj.updateDT(obj.PointStorage);
            end

            %% Perform current run

            while (obj.PointNum < obj.PointNumMax) || (obj.PointNumMax == 0)
                % Identify candidate simplexes
                pair = edges(obj.DT);
                cand_pair = pair(abs(obj.argDiff(pair(:,1), pair(:,2))) >= 2 * pi / 3,:);
                cand_simplex_ind = edgeAttachments(obj.DT, cand_pair);
                cand_simplex_ind = unique([cand_simplex_ind{:}]');
                cand_simplex = obj.DT(cand_simplex_ind,:);
                
                % Identify noncandidate simplexes
                noncand_simplex = obj.DT(:,:);
                noncand_simplex(cand_simplex_ind,:) = [];
                noncand_simplex = sort(noncand_simplex, 2);
                [prev_ind, prev_loc] = ismember(obj.SimplexGrad(:,1:3), noncand_simplex, 'rows');
                obj.SimplexGrad = obj.SimplexGrad(prev_ind,:);
                noncand_simplex(prev_loc(prev_ind),:) = [];
    
                % Calculate noncandidate simplexes absolute value gradient
                simplex_vector = permute(reshape(obj.DT.Points(reshape(noncand_simplex(:,2:end)', [], 1),:)', 2, 2, []), [3 1 2]) - ...
                        obj.DT.Points(noncand_simplex(:,1),:);
                grad_vector = permute(obj.absDiff(noncand_simplex(:,1), noncand_simplex(:,2:end)), [1 3 2]);
    
                grad_vector_x = simplex_vector;
                grad_vector_y = simplex_vector;
                grad_vector_x(:,1,:) = grad_vector;
                grad_vector_y(:,2,:) = grad_vector;
    
                obj.SimplexGrad = ...
                    [
                        obj.SimplexGrad
                        noncand_simplex ...
                        zeros([size(noncand_simplex,1) 1]) ...
                        [
                            obj.detND(grad_vector_x) ...
                            obj.detND(grad_vector_y)
                        ] ./ obj.detND(simplex_vector)
                    ];

                % Identify noncandidate simplex pairs
                noncand_ld_simplex = zeros([3 * size(obj.SimplexGrad, 1) 2]);
                permute_ind = nchoosek(1:3, 2);

                for i = 1:3
                    noncand_ld_simplex((1:size(obj.SimplexGrad, 1))+(i-1)*size(obj.SimplexGrad, 1),:) = obj.SimplexGrad(:,permute_ind(i,:));
                end

                [~, first_noncand_ld_simplex_ind, unique_noncand_ld_simplex_loc] = unique(noncand_ld_simplex, 'rows');
                second_noncand_ld_simplex_ind = setdiff(1:size(noncand_ld_simplex, 1), first_noncand_ld_simplex_ind)';
                first_noncand_ld_simplex_ind = first_noncand_ld_simplex_ind(unique_noncand_ld_simplex_loc(second_noncand_ld_simplex_ind));
                noncand_ld_simplex = noncand_ld_simplex(second_noncand_ld_simplex_ind,:);
                first_noncand_simplex_ind = mod(first_noncand_ld_simplex_ind-1, size(obj.SimplexGrad, 1)) + 1;
                second_noncand_simplex_ind = mod(second_noncand_ld_simplex_ind-1, size(obj.SimplexGrad, 1)) + 1;

                % Calculate noncandidate simplex pairs gradient flow
                noncand_ld_simplex_vector = obj.DT.Points(noncand_ld_simplex(:,2),:) - obj.DT.Points(noncand_ld_simplex(:,1),:);
                flow = ...
                    abs( ...
                        (obj.SimplexGrad(first_noncand_simplex_ind,end-1) - obj.SimplexGrad(second_noncand_simplex_ind,end-1)) .* obj.detND(noncand_ld_simplex_vector(:,1:end~=1,:)) - ...
                        (obj.SimplexGrad(first_noncand_simplex_ind,end) - obj.SimplexGrad(second_noncand_simplex_ind,end)) .* obj.detND(noncand_ld_simplex_vector(:,1:end~=2,:)) ...
                    );

                % Apply gradient flow upper bound
                flow_ind = (flow >= obj.FlowMax);

                % Weaken gradient flow upper bound if there are not enough triangulation points
                if (sum(flow_ind) == 0) && (obj.PointNum < obj.PointNumMin)
                    flow_ind = (flow >= max(flow) / 2);
                end

                % Identify gradient flow simplexes
                flow_simplex = ...
                    obj.SimplexGrad( ...
                        unique( ...
                            [
                                first_noncand_simplex_ind(flow_ind)
                                second_noncand_simplex_ind(flow_ind)
                            ] ...
                        ) ...
                    , 1:3);

                % Identify nonboundary simplexes
                nonboundary_simplex =...
                    [
                        cand_simplex
                        flow_simplex
                    ];

                % Identify nonboundary pairs
                permute_ind = nchoosek(1:3, 2);

                nonboundary_pair = nonboundary_simplex(:,permute_ind(1,:));
                nonboundary_pair_length = obj.coordDiff(nonboundary_simplex(:,permute_ind(1,1)), nonboundary_simplex(:,permute_ind(1,2)));

                for i = 2:3
                    next_nonboundary_pair_length = obj.coordDiff(nonboundary_simplex(:,permute_ind(i,1)), nonboundary_simplex(:,permute_ind(i,2)));
                    length_ind = (next_nonboundary_pair_length > nonboundary_pair_length);
                    nonboundary_pair(length_ind,:) = nonboundary_simplex(length_ind,permute_ind(i,:));
                    nonboundary_pair_length(length_ind) = next_nonboundary_pair_length(length_ind);
                end

                % Identify boundary simplexes
                boundary_simplex_ind = edgeAttachments(obj.DT, nonboundary_pair);
                boundary_simplex_ind = unique([boundary_simplex_ind{:}]');
                boundary_simplex = obj.DT(boundary_simplex_ind,:);
                boundary_simplex = boundary_simplex(sum(ismember(boundary_simplex, freeBoundary(obj.DT)), 2) == 2,:);

                % Identify boundary pairs
                permute_ind = nchoosek(1:3, 2);

                boundary_pair = boundary_simplex(:,permute_ind(1,:));
                boundary_pair_length = obj.coordDiff(boundary_simplex(:,permute_ind(1,1)), boundary_simplex(:,permute_ind(1,2)));

                for i = 2:3
                    next_boundary_pair_length = obj.coordDiff(boundary_simplex(:,permute_ind(i,1)), boundary_simplex(:,permute_ind(i,2)));
                    length_ind = (next_boundary_pair_length > boundary_pair_length);
                    boundary_pair(length_ind,:) = boundary_simplex(length_ind,permute_ind(i,:));
                    boundary_pair_length(length_ind) = next_boundary_pair_length(length_ind);
                end

                % Identify refinement simplexes
                ref_simplex =...
                    [
                        nonboundary_simplex
                        boundary_simplex
                    ];

                % Identify refinement pairs
                ref_pair =...
                    [
                        nonboundary_pair
                        boundary_pair
                    ];

                % Apply approximation tolerance condition
                ref_tol_ind = any(abs(obj.DT.Points(ref_pair(:,1),:) - obj.DT.Points(ref_pair(:,2),:)) > obj.Tol, 2);
                ref_pair = ref_pair(ref_tol_ind,:);

                simplex_tol_ind = ismember(obj.SimplexGrad(:,1:3), ref_simplex(~ref_tol_ind,:), 'rows') & ~obj.SimplexGrad(:,end-2);
                obj.SimplexGrad(simplex_tol_ind,end-2) = 1;

                % Identify refinement points
                ref_point = unique((obj.DT.Points(ref_pair(:,1),:) + obj.DT.Points(ref_pair(:,2),:)) / 2, 'rows');

                % Check whether there are no triangualtion points to add
                if isempty(ref_point) && ~any(simplex_tol_ind)
                    break
                end

                % Add new triangulation points
                obj.updateDT(ref_point);
            end

            %% Display algorithm progress

            if obj.Display == "on"
                fprintf( ...
                        [
                            '<strong>Refinement is completed.</strong>\n' ...
                            'Refinement time: <strong>%.2f</strong>\n\n'
                        ], ...
                        toc(obj.StartTime) ...
                    );

                obj.StartTime = tic;
            end

            %% Identify candidate regions and candidate points

            % Identify candidate pairs
            pair = edges(obj.DT);
            cand_pair = pair(abs(obj.argDiff(pair(:,1), pair(:,2))) >= 2 * pi / 3,:);

            % Set candidate regions dummy
            obj.CandRegion = {};

            % Set candidate points dummy
            obj.CandPoint = zeros(0,2);

            if isempty(cand_pair)
                warning('Neither zeros nor poles were found in the domain.');
            else
                % Identify candidate region pairs
                cand_simplex_ind = edgeAttachments(obj.DT, cand_pair);
                cand_simplex_ind = unique([cand_simplex_ind{:}]');

                permute_ind = nchoosek(1:3,2);
                permute_ind(2,:) = fliplr(permute_ind(2,:));

                ld_cand_simplex = zeros(3*length(cand_simplex_ind), 2);

                for i = 1:3
                    ld_cand_simplex((1:length(cand_simplex_ind))+(i-1)*length(cand_simplex_ind),:) = obj.DT(cand_simplex_ind,permute_ind(i,:));
                end

                cand_region_pair = ld_cand_simplex(~ismember(ld_cand_simplex, fliplr(ld_cand_simplex), 'rows'),:);

                % Set initial candidate region
                current_region = 1;
                current_pair = cand_region_pair(1,:);
                obj.CandRegion{current_region} = current_pair(1);
                cand_region_pair(1,:) = [];

                while ~isempty(cand_region_pair)
                    % Identify next candidate region pair
                    next_ind = 1:size(cand_region_pair, 1);
                    next_ind = next_ind(cand_region_pair(:,1) == current_pair(2));

                    % Check whether there are multiple next candidate region pairs
                    if sum(next_ind) > 1
                        % Identify the most outter next candidate region pair
                        current_vector = obj.DT.Points(current_pair(2),:) - obj.DT.Points(current_pair(1),:);
                        next_vector = obj.DT.Points(cand_region_pair(next_ind,2),:) - obj.DT.Points(cand_region_pair(next_ind,1),:);
                        [~, min_angle_ind] = ...
                            min(...
                                atan2(...
                                    current_vector(1) * next_vector(:,2) - current_vector(2) * next_vector(:,1),...
                                    current_vector(1) * next_vector(:,1) + current_vector(2) * next_vector(:,2)...
                                )...
                            );
                        next_ind = next_ind(min_angle_ind);
                    end

                    % Check whether there is next candidate region pair
                    if isempty(next_ind)
                        % Set next candidate region
                        current_region = current_region + 1;
                        current_pair = cand_region_pair(1,:);
                        obj.CandRegion{current_region} = current_pair(1);
                        cand_region_pair(1,:) = [];
                    else
                        % Add next candidate region pair to candidate region
                        current_pair = cand_region_pair(next_ind,:);
                        obj.CandRegion{current_region} =...
                            [
                                obj.CandRegion{current_region}
                                current_pair(1)
                            ];
                        cand_region_pair(next_ind,:) = [];
                    end
                end

                % Calculate argument change for each candidate region
                for i = 1:length(obj.CandRegion)
                    obj.CandPoint(i,1) = (obj.Domain(1,:) + obj.DomainNorm .* mean(obj.DT.Points(obj.CandRegion{i},:),1)) * [1; 1i];
                    arg_diff = obj.argDiff(obj.CandRegion{i}, circshift(obj.CandRegion{i},-1));
                    obj.CandPoint(i,2) = sum(arg_diff .* (arg_diff < 2 * pi / 3)) / (2 * pi);
                end      

                obj.CandPoint(:,2) = round(obj.CandPoint(:,2), ceil(-log10(sqrt(eps))));

                %% Displaying algorithm progress

                if obj.Display == "on"
                    fprintf( ...
                            [
                                '<strong>Classification is completed.</strong>\n' ...
                                'Classification time: <strong>%0.2f</strong> s\n\n'
                            ], ...
                            toc(obj.StartTime) ...
                        );

                    disp(sortrows(array2table(obj.CandPoint, VariableNames=["z" "k"]), 'k', 'descend'));
                end
            end
        end

        function visTriang(obj, region_type)
            %% visTriang: triangulation visualization
            %
            % INPUT
            %
            %   (optional) region_type: region type,
            %                           subset of [0 -1 1]
            %   default: [-1 1]
            %
            % OUTPUT
            %   

            arguments
                obj
                region_type (1,:) {mustBeMember(region_type, [0 -1 1])} = [-1 1]
            end

            % Set color set
            color_set = ...
                [
                    0.8500 0.3250 0.0980;
                    0 0 0;
                    0 0.4470 0.7410
                ];
    
            %% Plot triangulation

            hold on

            % Construct renormalized triangulation
            TR = triangulation(obj.DT(:,:), obj.Domain(1,:) + obj.DomainNorm .* obj.DT.Points);

            % Plot triangulation
            triplot( ...
                    TR, ...
                    '-k' ...
                );

            % Plot candidate regions and candidate points
            for type = region_type
                % Identify candidate points with a certain type
                if type == 0
                    ind = (obj.CandPoint(:,end) == 0);
                else
                    ind = (type * obj.CandPoint(:,end) > 0);
                end

                % Plot candidate points
                scatter(...
                        real(obj.CandPoint(ind,1)),...
                        imag(obj.CandPoint(ind,1)),...
                        [],...
                        color_set(type+2,:),...
                        'filled' ...
                    )

                % Set candidate region
                cand_region = obj.CandRegion(ind);

                % Plot candidate region
                for i = 1:length(cand_region)
                    cand_region_point = obj.Domain(1,:) + obj.DomainNorm .* obj.DT.Points(cand_region{i},:);
                    cand_region_point =...
                        [
                            cand_region_point
                            cand_region_point(1,:)
                        ];

                    plot(...
                            cand_region_point(:,1), ...
                            cand_region_point(:,2), ...
                            Color=color_set(type+2,:), ...
                            LineWidth=get(0, 'defaultLineLineWidth') + 1 ...
                        )
                end
            end

            hold off

            % Set axes limits
            axis(obj.Domain(:));
    
            % Set x-label
            xlabel( ...
                    '$x$', ...
                    Interpreter='latex' ...
                );

            % Set y-label
            ylabel( ...
                    '$y$', ...
                    Interpreter='latex' ...
                );
        end
    end

    methods (Access=private, Hidden)
        function updateDT(obj, ref_point)
            %% updateDT: Delaunay triangulation update
            %
            % INPUT
            %
            %   ref_point: refinement points,
            %              ?x2 matrix
            %
            % OUTPUT
            %   

            % Apply maximum number of triangulation points condition
            if obj.PointNumMax > 0
                point_max = min(size(ref_point,1), obj.PointNumMax - obj.PointNum);
                obj.PointStorage = ref_point(point_max+1:size(ref_point, 1),:);
                ref_point = ref_point(1:point_max,:);
            end

            % Add new triangulation points
            obj.DT.Points =...
                [
                    obj.DT.Points
                    ref_point
                ];

            % Set number of triangulation points
            obj.PointNum = size(obj.DT.Points, 1);

            % Check whether batch size is nonzero
            if obj.BatchSize == 0
                % Calculate new function values
                obj.FuncEval = ...
                    [
                        obj.FuncEval
                        obj.Func(ref_point)
                    ];
            else
                % Calculate batched function values
                ref_func_eval = zeros(size(ref_point,1),1);

                for i = 1:obj.BatchSize:length(ref_func_eval)
                    batch_ind = i:min(i+obj.BatchSize-1,length(ref_func_eval));
                    ref_func_eval(batch_ind) = obj.Func(ref_point(batch_ind,:));
                end

                % Add new function values
                obj.FuncEval =...
                    [
                        obj.FuncEval
                        ref_func_eval
                    ];
            end

            %% Displaying algorithm progress
    
            if obj.Display == "on"
                fprintf( ...
                        [
                            'Points number: <strong>%i</strong>\n' ...
                            'Refinement time: <strong>%.2f</strong> s\n\n'
                        ], ...
                        obj.PointNum, ...
                        toc(obj.StartTime) ...
                    );
            end
        end

        function value = absDiff(obj, first_point, second_point)
            %% absDiff: absolute value difference in a pair
            %
            % INPUT
            %
            %   first_point: first point,
            %                vector
            %
            %   second_point: second point,
            %                 vector
            %
            % OUTPUT
            %
            %   value: absolute value difference,
            %          vector

            value = log(abs(obj.FuncEval(second_point) ./ obj.FuncEval(first_point)));
        end

        function value = argDiff(obj, first_point, second_point)
            %% argDiff: argument difference in a pair
            %
            % INPUT
            %
            %   first_point: first point,
            %                vector
            %
            %   second_point: second point,
            %                 vector
            %
            % OUTPUT
            %
            %   value: argument difference,
            %          vector

            arg = obj.FuncEval(second_point) ./ obj.FuncEval(first_point);
            value = angle(arg);
            value((arg == 0) | isinf(arg) | isnan(arg)) = 2 * pi / 3;
        end

        function value = coordDiff(obj, first_point, second_point)
            %% coordDiff: coordinate difference in a pair
            %
            % INPUT
            %
            %   first_point: first point,
            %                vector
            %
            %   second_point: second point,
            %                 vector
            %
            % OUTPUT
            %
            %   value: coordinate difference,
            %          vector

            value = vecnorm(obj.DT.Points(second_point,:) - obj.DT.Points(first_point,:), 2, 2);
        end

        function value = detND(obj, M)
            %% detND: vectorized matrix determinant
            %
            % INPUT
            %
            %   M: 3-dimensional array
            %
            % OUTPUT
            %
            %   value: determinant,
            %          vector

            if prod(size(M, [2 3])) <= 1
                value = M;
            else
                value = 0;
                
                for i = 1:size(M, 2)
                    value = value + (-1)^(i+1) * M(:,i,1) .* obj.detND(M(:,1:size(M,2) ~= i,2:end));
                end
            end
        end
    end
end