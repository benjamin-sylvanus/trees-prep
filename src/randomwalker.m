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
            obj.step = step;
%             obj = obj.initPosition(flag, bounds(1), bounds(2), bounds(3), ...
%             swc, LUT, A, pairs);
            obj.iter = iter;
            obj.chunkSize = chunkSize;
            obj.flag = false;
%             obj = obj.init_rands();
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

        function inside = checkpos(obj, pos, swc, LUT, A, pairs)

            % check lookup-table using subscripted index of [x y z]
            ps = ceil(pos);
            indicies = LUT(ps(1), ps(2), ps(3));
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
                        inside = obj.pointbasedswc2v(x0, y0, z0, ...
                            x1, x2, y1, ...
                            y2, z1, z2, ...
                            r1, r2, obj.step, dist) | inside;
                    else
                        inside = obj.pointbasedswc2v(x0, y0, z0, ...
                            x2, x1, y2, ...
                            y1, z2, z1, ...
                            r2, r1, obj.step, dist) | inside;
                    end

                end

            end

        end

    end

    methods (Static)

        function pos = randxyz(sx, sy, sz)
            %             x = randi([round(sx / 4), round(3 * sx / 4)], 1000, 1); y = randi([round(sy / 4), round(3 * sy / 4)], 1000, 1); z = randi([round(sz / 4), round(3 * sz / 4)], 1000, 1);
            x = randi(sx-1, 10000, 1) + 2 .* rand(10000,1) - 1; y = randi(sy-1, 10000, 1) + 2 .* rand(10000,1) - 1; z = randi(sz-1, 10000, 1) + 2 .* rand(10000,1) - 1;
            pos = [x, y, z];
        end

        function pos = pointbasedswc2v(x0, y0, z0, x1, x2, y1, y2, z1, z2, r1, r2, step ,dist)


        t = ((x0 - x1) * (x2 - x1) + (y0 - y1) * (y2 - y1) + (z0 - z1) * (z2 - z1)) ./ ...
            (dist);
        x = x1 + (x2 - x1) * t;
        y = y1 + (y2 - y1) * t;
        z = z1 + (z2 - z1) * t;

%         if (x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2 < r1^2
        if dist < r1^2
            list1 = false;
        else
            list1 = (x - x1) .* (x - x2) + (y - y1) .* (y - y2) + (z - z1) .* (z - z2) < 0;
        end

        if list1
            dist2 = (x0 - x) .^ 2 + (y0 - y) .^ 2 + (z0 - z) .^ 2;

            %     r = r1 + sqrt((x-x1).^2 + (y-y1).^2 + (z-z1).^2) /...
            %         sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) * (r2-r1);

            %     r = ( c + r2 ) / (sqrt ( 1 - ( |r1-r2 | / l ) )

            %     c = ( |r1 - r2| * l ) / L
            rd = abs(r1 - r2);

            % distance from orthogonal vector to p2
            l = sqrt((x - x2) .^ 2 + (y - y2) .^ 2 + (z - z2) .^ 2);

            % distance from p1 -> p2
            L = sqrt(dist);

            c = (rd * l) ./ L;
            r = (c + r2) ./ sqrt(1 - ((rd / L) .^ 2));
            pos1 = dist2 < ((r-step) .^ 2); % smaller in one line and less than and equal
            pos = pos1;
        else
            pos2 = (((x0 - x1) .^ 2 + (y0 - y1) .^ 2 + (z0 - z1) .^ 2) < ((r1-step) ^ 2)) | ...
                (((x0 - x2) .^ 2 + (y0 - y2) .^ 2 + (z0 - z2) .^ 2) < ((r2-step) ^ 2));
            pos = pos2;
        end

end

    end

end
