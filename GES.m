classdef GES < handle
    properties (SetAccess=private)
        %% Domain: analyzed domain
        Domain (2,2)
        %   matrix

        %% DomainNorm: domain normalization
        DomainNorm (1,1)
        %   scalar

        %% Tol: approximation relative tolerance level
        Tol (1,1)
        %   scalar

        %% Func: analyzed function
        Func function_handle
        %   function handle

        %% PointNumMin: triangulation points number minimum
        PointNumMin (1,1)
        %   scalar

        %% PropMax: absolute value flow maximum
        PropMax (1,1)
        %   scalar

        %% DT: Delaunay triangulation
        DT = delaunayTriangulation()
        %   Delaunay triangulation class

        %% FuncEval(i): function evaluations
        FuncEval (:,1)
        %   vector

        %% PointNum: triangulation points number
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
        %% PointNumMax: triangulation points number maximum
        PointNumMax (1,1) {mustBeInteger, mustBeGreaterThan(PointNumMax, 2)} = 3
        %   scalar

        %% BatchSize: batch size
        BatchSize (1,1) {mustBeInteger, mustBeNonnegative} = 0
        %   scalar

        %% Display: output flag
        Display {mustBeMember(Display, ["off" "on"])} = 'off'
        %   either "off" or "on"
    end

    properties (Access=private, Hidden)
        %% SimplexProp(i,j): simplex logarithmic properties
        SimplexProp (:,:) = zeros(0,6)
        %   matrix

        %% PointStorage(i,j): leftover points storage
        PointStorage (:,:)
        %   matrix

        %% StartTime: start time
        StartTime (1,1)
        %   scalar
    end

    methods
        function obj = GES(Func, Domain, Tol, PointNumMax, options)
            %% Global Equation Solver (GES)
            %   Func - the analyzed function, function handle f(z)
            %   Domain(i,j) - the analyzed domain boundaries, a matrix 2x2
            %   Tol - approximation relative tolerance level, a
            %         positive scalar with value less than 1       
            %   PointNumMax - triangulation points number maximum, an 
            %                 integer greater than 2
            %   PointNumMin (optional) - triangulation points number
            %                            minimum, a positive integer
            %                            (default 25)
            %   PropMax (optional) - absolute value flow maximum, a 
            %                        positive scalar (default 1)
            %   BatchSize (optional) - batch size, a nonnegative integer
            %                          (default 0, no batching)
            %   Display (optional) - output flag, either "off" or "on" 
            %                        (default "off")
        
            arguments
                Func function_handle
                Domain (2,2) {mustBeReal, mustBeDimSorted(Domain, 1), mustBeNonzeroDifference}
                Tol (1,1) {mustBePositive, mustBeLessThan(Tol, 1)}
                PointNumMax (1,1) {mustBeInteger, mustBeGreaterThan(PointNumMax, 2)}
                options.PointNumMin (1,1) {mustBeInteger, mustBePositive} = 25
                options.PropMax (1,1) {mustBePositive} = 1
                options.BatchSize (1,1) {mustBeInteger, mustBeNonnegative} = 0
                options.Display {mustBeMember(options.Display, ["off" "on"])} = 'off'
            end

            %% Setting class properties

            obj.Domain = Domain;

            side = diff(Domain,1,1);
            obj.DomainNorm = min(side);
            side_rel = side / obj.DomainNorm;

            n_x = ceil(side_rel(1)) + 1;
            n_y = ceil(side_rel(2)) + 1;
            [initial_point_x, initial_point_y] = ndgrid(...
                        linspace(0, side_rel(1), n_x),...
                        linspace(0, side_rel(2), n_y)...
                    );
            initial_point =...
                [
                    initial_point_x(:) initial_point_y(:)
                ];
            initial_point_x = initial_point_x(1:end-1,1:end-1);
            initial_point_y = initial_point_y(1:end-1,1:end-1);
            initial_point =...
                [
                    initial_point
                    [initial_point_x(:) initial_point_y(:)] + side_rel ./ [n_x - 1 n_y - 1] / 2
                ];

            obj.PointStorage = initial_point;

            obj.Tol = Tol;

            obj.Func = @(z) Func((obj.Domain(1,:) + obj.DomainNorm .* z) * [1; 1i]);

            obj.PointNumMax = PointNumMax;
            obj.PointNumMin = options.PointNumMin;
            obj.PropMax = options.PropMax;

            obj.BatchSize = options.BatchSize;
            obj.Display = options.Display;

            obj.StartTime = tic;

            obj.fitTriang;
        end

        function fitTriang(obj)
            %% fitTriang: triangulation fitting

            %% Completing previous run

            if ~isempty(obj.PointStorage)
                obj.updateDT(obj.PointStorage);
            end

            %% Performing current run

            while obj.PointNum < obj.PointNumMax
                %% Identifying candidate simplexes

                pair = edges(obj.DT);
                cand_pair = pair(abs(obj.argDiff(pair(:,1), pair(:,2))) >= 2 * pi / 3,:);
                cand_simplex_ind = edgeAttachments(obj.DT, cand_pair);
                cand_simplex_ind = unique([cand_simplex_ind{:}]');
                cand_simplex = obj.DT(cand_simplex_ind,:);
                
                %% Updating properties of noncandidate simplexes
    
                noncand_simplex = obj.DT(:,:);
                noncand_simplex(cand_simplex_ind,:) = [];
                noncand_simplex = sort(noncand_simplex,2);
                [prev_ind, prev_loc] = ismember(obj.SimplexProp(:,1:3), noncand_simplex, 'rows');
                obj.SimplexProp = obj.SimplexProp(prev_ind,:);
                noncand_simplex(prev_loc(prev_ind),:) = [];
    
                simplex_vector = permute(reshape(obj.DT.Points(reshape(noncand_simplex(:,2:end)',[],1),:)',2,2,[]),[3 1 2]) -...
                        obj.DT.Points(noncand_simplex(:,1),:);
                prop_vector = permute(obj.absDiff(noncand_simplex(:,1), noncand_simplex(:,2:end)),[1 3 2]);
    
                prop_vector_x = simplex_vector;
                prop_vector_y = simplex_vector;
                prop_vector_x(:,1,:) = prop_vector;
                prop_vector_y(:,2,:) = prop_vector;
    
                obj.SimplexProp =...
                    [
                        obj.SimplexProp
                        noncand_simplex...
                        zeros(size(noncand_simplex,1),1)...
                        [
                            obj.detND(prop_vector_x)...
                            obj.detND(prop_vector_y)
                        ] ./ obj.detND(simplex_vector)
                    ];

                %% Identifying noncandidate simplex pairs

                noncand_ld_simplex = zeros(3*size(obj.SimplexProp,1),2);
                permute_ind = nchoosek(1:3,2);

                for i = 1:3
                    noncand_ld_simplex((1:size(obj.SimplexProp,1))+(i-1)*size(obj.SimplexProp,1),:) = obj.SimplexProp(:,permute_ind(i,:));
                end

                [~, first_noncand_ld_simplex_ind, unique_noncand_ld_simplex_loc] = unique(noncand_ld_simplex, 'rows');
                second_noncand_ld_simplex_ind = setdiff(1:size(noncand_ld_simplex,1),first_noncand_ld_simplex_ind)';
                first_noncand_ld_simplex_ind = first_noncand_ld_simplex_ind(unique_noncand_ld_simplex_loc(second_noncand_ld_simplex_ind));
                noncand_ld_simplex = noncand_ld_simplex(second_noncand_ld_simplex_ind,:);
                first_noncand_simplex_ind = mod(first_noncand_ld_simplex_ind-1,size(obj.SimplexProp,1)) + 1;
                second_noncand_simplex_ind = mod(second_noncand_ld_simplex_ind-1,size(obj.SimplexProp,1)) + 1;

                %% Ranking noncandidate simplex pairs

                noncand_ld_simplex_vector = obj.DT.Points(noncand_ld_simplex(:,2),:) - obj.DT.Points(noncand_ld_simplex(:,1),:);
                prop =...
                    abs(...
                        (obj.SimplexProp(first_noncand_simplex_ind,end-1) - obj.SimplexProp(second_noncand_simplex_ind,end-1)) .* obj.detND(noncand_ld_simplex_vector(:,1:end~=1,:)) -...
                        (obj.SimplexProp(first_noncand_simplex_ind,end) - obj.SimplexProp(second_noncand_simplex_ind,end)) .* obj.detND(noncand_ld_simplex_vector(:,1:end~=2,:))...
                    );

                prop_ind = (prop >= obj.PropMax);

                if (sum(prop_ind) == 0) && (obj.PointNum < obj.PointNumMin)
                    prop_ind = (prop >= max(prop) / 2);
                end

                prop_simplex = obj.SimplexProp(unique(...
                    [
                        first_noncand_simplex_ind(prop_ind)
                        second_noncand_simplex_ind(prop_ind)
                    ]),1:3);

                %% Identifying nonboundary pairs

                nonboundry_simplex =...
                    [
                        cand_simplex
                        prop_simplex
                    ];

                permute_ind = nchoosek(1:3,2);

                nonboundary_pair = nonboundry_simplex(:,permute_ind(1,:));
                nonboundary_pair_length = obj.coordDiff(nonboundry_simplex(:,permute_ind(1,1)), nonboundry_simplex(:,permute_ind(1,2)));

                for i = 2:3
                    next_nonboundary_pair_length = obj.coordDiff(nonboundry_simplex(:,permute_ind(i,1)), nonboundry_simplex(:,permute_ind(i,2)));
                    length_ind = (next_nonboundary_pair_length > nonboundary_pair_length);
                    nonboundary_pair(length_ind,:) = nonboundry_simplex(length_ind,permute_ind(i,:));
                    nonboundary_pair_length(length_ind) = next_nonboundary_pair_length(length_ind);
                end

                %% Identifying boundary pairs

                boundary_simplex_ind = edgeAttachments(obj.DT, nonboundary_pair);
                boundary_simplex_ind = unique([boundary_simplex_ind{:}]');
                boundary_simplex = obj.DT(boundary_simplex_ind,:);
                boundary_simplex = boundary_simplex(sum(ismember(boundary_simplex, freeBoundary(obj.DT)),2) == 2,:);

                permute_ind = nchoosek(1:3,2);

                boundary_pair = boundary_simplex(:,permute_ind(1,:));
                boundary_pair_length = obj.coordDiff(boundary_simplex(:,permute_ind(1,1)), boundary_simplex(:,permute_ind(1,2)));

                for i = 2:3
                    next_boundary_pair_length = obj.coordDiff(boundary_simplex(:,permute_ind(i,1)), boundary_simplex(:,permute_ind(i,2)));
                    length_ind = (next_boundary_pair_length > boundary_pair_length);
                    boundary_pair(length_ind,:) = boundary_simplex(length_ind,permute_ind(i,:));
                    boundary_pair_length(length_ind) = next_boundary_pair_length(length_ind);
                end

                %% Verifying tolerance condition

                ref_simplex =...
                    [
                        nonboundry_simplex
                        boundary_simplex
                    ];

                ref_pair =...
                    [
                        nonboundary_pair
                        boundary_pair
                    ];

                ref_tol_ind = any(abs(obj.DT.Points(ref_pair(:,1),:) - obj.DT.Points(ref_pair(:,2),:)) > obj.Tol,2);
                ref_pair = ref_pair(ref_tol_ind,:);

                simplex_tol_ind = ismember(obj.SimplexProp(:,1:3), ref_simplex(~ref_tol_ind,:), 'rows') & ~obj.SimplexProp(:,end-2);
                obj.SimplexProp(simplex_tol_ind,end-2) = 1;

                %% Adding new triangulation points

                ref_point = unique((obj.DT.Points(ref_pair(:,1),:) + obj.DT.Points(ref_pair(:,2),:)) / 2, 'rows');

                if isempty(ref_point) && ~any(simplex_tol_ind)
                    break
                end

                obj.updateDT(ref_point);
            end

            %% Displaying algorithm progress

            if obj.Display == "on"
                char_num = fprintf('<strong>Triangulation is completed.</strong>\n');
                fprintf('Elapsed time: <strong>%0.2f</strong> s\n', toc(obj.StartTime));
                fprintf([repmat('-',1,char_num) '\n']);
            end

            %% Identifying candidate pairs

            pair = edges(obj.DT);
            cand_pair = pair(abs(obj.argDiff(pair(:,1), pair(:,2))) >= 2 * pi / 3,:);

            %% Identifying candidate regions

            obj.CandRegion = {};
            obj.CandPoint = zeros(0,2);

            if isempty(cand_pair)
                warning('Neither zeros nor poles were found in the domain.');
            else
                %% Identifying candidate region boundaries

                cand_simplex_ind = edgeAttachments(obj.DT, cand_pair);
                cand_simplex_ind = unique([cand_simplex_ind{:}]');

                permute_ind = nchoosek(1:3,2);
                permute_ind(2,:) = fliplr(permute_ind(2,:));

                ld_cand_simplex = zeros(3*length(cand_simplex_ind),2);

                for i = 1:3
                    ld_cand_simplex((1:length(cand_simplex_ind))+(i-1)*length(cand_simplex_ind),:) = obj.DT(cand_simplex_ind,permute_ind(i,:));
                end

                cand_region_pair = ld_cand_simplex(~ismember(ld_cand_simplex, fliplr(ld_cand_simplex), 'rows'),:);

                %% Classifying candidate regions

                current_region = 1;
                current_pair = cand_region_pair(1,:);
                obj.CandRegion{current_region} = current_pair(1);
                cand_region_pair(1,:) = [];

                while ~isempty(cand_region_pair)
                        next_ind = 1:size(cand_region_pair,1);
                        next_ind = next_ind(cand_region_pair(:,1) == current_pair(2));

                        if sum(next_ind) > 1
                            current_vector = obj.DT.Points(current_pair(2),:) - obj.DT.Points(current_pair(1),:);
                            next_vector = obj.DT.Points(cand_region_pair(next_ind,2),:) - obj.DT.Points(cand_region_pair(next_ind,1),:);
                            [~, min_angle_ind] = min(...
                                    atan2(...
                                        current_vector(1) * next_vector(:,2) - current_vector(2) * next_vector(:,1),...
                                        current_vector(1) * next_vector(:,1) + current_vector(2) * next_vector(:,2)...
                                    )...
                                );
                            next_ind = next_ind(min_angle_ind);
                        end

                        if isempty(next_ind)
                            current_region = current_region + 1;

                            current_pair = cand_region_pair(1,:);
                            obj.CandRegion{current_region} = current_pair(1);
                            cand_region_pair(1,:) = [];
                        else
                            current_pair = cand_region_pair(next_ind,:);

                            obj.CandRegion{current_region} =...
                                [
                                    obj.CandRegion{current_region}
                                    current_pair(1)
                                ];
                            cand_region_pair(next_ind,:) = [];
                        end
                end

                for i = 1:length(obj.CandRegion)
                    obj.CandPoint(i,1) = (obj.Domain(1,:) + obj.DomainNorm .* mean(obj.DT.Points(obj.CandRegion{i},:),1)) * [1; 1i];
                    obj.CandPoint(i,2) = sum(obj.argDiff(obj.CandRegion{i}, circshift(obj.CandRegion{i},-1))) / (2 * pi);
                end      

                obj.CandPoint(:,2) = round(obj.CandPoint(:,2), ceil(-log10(sqrt(eps))));
            end
        end

        function visTriang(obj, RegionType)
            %% visTriang: triangulation visualization
            %   RegionType (optional) - region type, a member of [0 -1 1]

            arguments
                obj
                RegionType (1,:) {mustBeMember(RegionType, [0 -1 1])} = [0 -1 1]
            end

            %% Setting figure parameters

            set_color =...
                    [
                    0.8500 0.3250 0.0980;
                    0 0 0;
                    0 0.4470 0.7410
                    ];
    
            figure(Units='normalized', Position=[0.2 0.2 0.6 0.6]);

            %% Plotting triangulation

            hold on

            TR = triangulation(obj.DT(:,:), obj.Domain(1,:) + obj.DomainNorm .* obj.DT.Points);
            triplot(TR, '-k');

            for type = RegionType
                if type == 0
                    ind = (obj.CandPoint(:,end) == 0);
                else
                    ind = (type * obj.CandPoint(:,end) > 0);
                end

                scatter(...
                        real(obj.CandPoint(ind,1)),...
                        imag(obj.CandPoint(ind,1)),...
                        25,...
                        set_color(type+2,:),...
                        'filled' ...
                    )

                cand_region = obj.CandRegion(ind);

                for i = 1:length(cand_region)
                    cand_region_point = obj.Domain(1,:) + obj.DomainNorm .* obj.DT.Points(cand_region{i},:);
                    cand_region_point =...
                        [
                            cand_region_point
                            cand_region_point(1,:)
                        ];

                    plot(...
                            cand_region_point(:,1),...
                            cand_region_point(:,2),...
                            Color=set_color(type+2,:),...
                            LineWidth=2 ...
                        )
                end
            end

            hold off

            axis(obj.Domain(:));
    
            set(gca, 'FontSize', 16);
    
            xlabel('$x$', Interpreter='latex', FontSize=18);
            ylabel('$y$', Interpreter='latex', FontSize=18);
        end
    end

    methods (Access=private, Hidden)
        function updateDT(obj, RefPoint)
            %% updateDT: Delaunay triangulation update
            %   NewPoint(i,j) - new triangulation points, matrix

            point_max = min(size(RefPoint,1), obj.PointNumMax - obj.PointNum);
            obj.PointStorage = RefPoint(point_max+1:size(RefPoint,1),:);
            ref_point = RefPoint(1:point_max,:);

            obj.DT.Points =...
                [
                    obj.DT.Points
                    ref_point
                ];
            obj.PointNum = size(obj.DT.Points,1);

            if obj.BatchSize == 0
                obj.FuncEval =...
                    [
                        obj.FuncEval
                        obj.Func(ref_point)
                    ];
            else
                ref_func_eval = zeros(size(ref_point,1),1);

                for i = 1:obj.BatchSize:length(ref_func_eval)
                    batch_ind = i:min(i+obj.BatchSize-1,length(ref_func_eval));
                    ref_func_eval(batch_ind) = obj.Func(ref_point(batch_ind,:));
                end

                obj.FuncEval =...
                    [
                        obj.FuncEval
                        ref_func_eval
                    ];
            end

            %% Displaying algorithm progress
    
            if obj.Display == "on"
                char_num = fprintf('Number of triangulation points: <strong>%i</strong>\n', obj.PointNum);
                fprintf('Elapsed time: <strong>%0.2f</strong> s\n', toc(obj.StartTime));
                fprintf([repmat('-',1,char_num) '\n']);
            end
        end

        function value = absDiff(obj, FirstPoint, SecondPoint)
            %% absDiff: absolute value difference in a pair
            %   FirstPoint(i) - first point, vector
            %   SecondPoint(i) - second point, vector

            value = log(abs(obj.FuncEval(SecondPoint) ./ obj.FuncEval(FirstPoint)));
        end

        function value = argDiff(obj, FirstPoint, SecondPoint)
            %% argDiff: argument difference in a pair
            %   FirstPoint(i) - first point, vector
            %   SecondPoint(i) - second point, vector

            arg = obj.FuncEval(SecondPoint) ./ obj.FuncEval(FirstPoint);
            value = angle(arg);
            value((arg == 0) | isinf(arg) | isnan(arg)) = 2 * pi / 3;
        end

        function value = coordDiff(obj, FirstPoint, SecondPoint)
            %% coordDiff: coordinate difference in a pair
            %   FirstPoint(i) - first point, vector
            %   SecondPoint(i) - second point, vector

            value = vecnorm(obj.DT.Points(SecondPoint,:) - obj.DT.Points(FirstPoint,:),2,2);
        end

        function value = detND(obj, M)
            %% detND: vectorized matrix determinant
            %   M(i,j,k) - array

            if prod(size(M,[2 3])) <= 1
                value = M;
            else
                value = 0;
                for i = 1:size(M,2)
                    value = value + (-1)^(i+1) * M(:,i,1) .* obj.detND(M(:,1:size(M,2) ~= i,2:end));
                end
            end
        end
    end
end

function mustBeDimSorted(x, dim)
    if ~issorted(x,dim)
        eid = 'Data:notSorted';
        msg = 'Rows of the argument must be sorted.';
        error(eid, msg);
    end
end

function mustBeNonzeroDifference(x)
    if any(diff(x,1,1) == 0)
        eid = 'Data:zeroDifference';
        msg = 'Argument must have different values in columns.';
        error(eid, msg);
    end
end