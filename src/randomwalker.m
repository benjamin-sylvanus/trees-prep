classdef randomwalker
    %RANDOMWALKER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        curr;
        next;
        step;
    end

    methods

        function obj = randomwalker(flag, step, bounds, swc, LUT, A, pairs)
            %RANDOMWALKER Construct an instance of this class
            %   Detailed explanation goes here

            obj = obj.initPosition(flag, bounds(1), bounds(2), bounds(3), ...
            swc, LUT, A, pairs);
            obj.step = step;
        end

        function obj = initPosition(obj, flag, sx, sy, sz, swc, LUT, A, pairs)

            if flag
                outside = true;

                while outside
                    % get a random position
                    pos = obj.randxyz(sx, sy, sz);
                    % check that position is inside the cell
                    inside = checkpos(obj, pos, swc, LUT, A, pairs);

                    if inside
                        outside = false;
                    end

                end

                obj.curr = pos;
            else
                obj.curr = randxyx(sx, sy, sz);
            end

        end

        function obj = setnext(obj)
            vect = random_unit_vector(3, 1);
            delta = vect * obj.step;
            pos2 = obj.curr + delta;
            obj.next = pos2;
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
                    p1 = swc{baseid, 2:5};
                    p2 = swc{targetid, 2:5};
                    x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
                    x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);
                    inside = pointbasedswc2v(x0, y0, z0, x1, x2, y1, y2, z1, z2, r1, r2) | inside;

                end

            end

        end

        function pos = randxyz(obj, sx, sy, sz)
            x = randi(sx, 1); y = randi(sy, 1); z = randi(sz, 1);
            pos = [x; y; z];
        end

    end

end
