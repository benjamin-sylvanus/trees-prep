classdef randomwalker
    %RANDOMWALKER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        curr;
        next;
        step;
        iter;
        rxyz;
    end

    methods (Access = public)

        function obj = randomwalker(flag, step, bounds, swc, LUT, A, pairs, iter)
            %RANDOMWALKER Construct an instance of this class
            %   Detailed explanation goes here

            obj = obj.initPosition(flag, bounds(1), bounds(2), bounds(3), ...
            swc, LUT, A, pairs);
            obj.step = step;
            obj.iter = iter;
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

                        % x0 = ps(1); y0 = ps(2); z0 = ps(3);

                        %                     bx = ~ismember(x0,42:55);
                        %                     by = ~ismember(y0,42:55);
                        %                     bz = ~ismember(z0,10:35);
                        %                     bcombined = bx |  by | bz;
                        %
                        %                     if inside(1) && bcombined
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
            obj.rxyz = random_unit_vector(3, obj.iter);

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
                    inside = pointbasedswc2v(x0, y0, z0, x1, x2, y1, y2, z1, z2, r1, r2, zeros(1, 1, "logical")) | inside;

                end

            end

        end

    end

    methods (Static)

        function pos = randxyz(sx, sy, sz)
            x = randi([round(sx / 4), round(3 * sx / 4)], 1000, 1); y = randi([round(sy / 4), round(3 * sy / 4)], 1000, 1); z = randi([round(sz / 4), round(3 * sz / 4)], 1000, 1);
            pos = [x, y, z];
        end

    end

end
