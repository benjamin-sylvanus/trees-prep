classdef randomwalker
    %RANDOMWALKER Summary of this class goes here
    %   Detailed explanation goes here

    properties

        curr; % current position [x;y;z]
        next; % next position [x;y;z]
        step; % step size
        iter; % ^ depreciated total_iterations in simulation
        rxyz; % N random xyz vectors -> 3 x chunkSize
        chunkSize; % sub-unit of iteration. store values for current chunk only
        flag; % on failed exit i.e. collision we stop for 1 step.
    end

    methods (Access = public)

        function obj = randomwalker(flag, step, bounds, swc, LUT, A, pairs, iter, chunkSize)
            %RANDOMWALKER Construct an instance of this class
            %   Detailed explanation goes here

            obj = obj.initPosition(flag, bounds(1), bounds(2), bounds(3), ...
            swc, LUT, A, pairs);
            obj.step = step;
            obj.iter = iter;
            obj.chunkSize = chunkSize;
            obj.flag = false;
            obj = obj.init_rands();
        end

        function obj = initPosition(obj, flag, sx, sy, sz, swc, LUT, A, pairs)

            if flag
                outside = true;

                while outside
                    % get a random position
                    pos = obj.randxyz(sx, sy, sz);

                    for i = 1:size(pos, 1)
                        ps = pos(i, :)';
                        % check that position is inside the cell
                        inside = checkpos(obj, ps, swc, LUT, A, pairs);

                        if inside(1)
                            pos = ps;
                            outside = false;
                            break;
                        end

                    end

                end

                obj.curr = pos;
            else
                obj.curr = randxyz(sx, sy, sz);
            end

        end

        function obj = init_rands(obj)
%             obj.rxyz = random_unit_vector(3, obj.iter);
        end

        function next = setnext(obj, i)
            vect = obj.rxyz(:, i);
            delta = vect * obj.step;
            next = obj.curr + delta;
        end

        function inside = checkpos(~, pos, swc, LUT, A, pairs)

            % check lookup-table using subscripted index of [x y z]
            indicies = LUT(pos(1), pos(2), pos(3));
            x0 = pos(1); y0 = pos(2); z0 = pos(3);
            inside = false;

            if (indicies ~= 0)

                % get pairs from A
                pairlist = A{indicies};
                ps = pairs(pairlist, :);

                % pairlist => [child, parent]
                children = ps(:, 1);
                parents = ps(:, 2);

                % for each pair: check if point is inside
                for i = 1:length(children)

                    % get base and target ids
                    baseid = children(i); targetid = parents(i);
                    p1 = swc(baseid, 2:5);
                    p2 = swc(targetid, 2:5);
                    x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
                    x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);
                    dist = (x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2;

                    if r1 > r2
                        inside = pointbasedswc2v(x0, y0, z0, ...
                            x1, x2, y1, ...
                            y2, z1, z2, ...
                            r1, r2, zeros(1, 1, "logical"), dist) | inside;
                    else
                        inside = pointbasedswc2v(x0, y0, z0, ...
                            x2, x1, y2, ...
                            y1, z2, z1, ...
                            r2, r1, zeros(1, 1, "logical"), dist) | inside;
                    end

                end

            end

        end

    end

    methods (Static)

        function pos = randxyz(sx, sy, sz)
%             x = randi([round(sx / 4), round(3 * sx / 4)], 1000, 1); y = randi([round(sy / 4), round(3 * sy / 4)], 1000, 1); z = randi([round(sz / 4), round(3 * sz / 4)], 1000, 1);
            x = randi(sx,10000,1); y = randi(sy,10000,1); z = randi(sz,10000,1);
            pos = [x, y, z];
        end

    end

end
