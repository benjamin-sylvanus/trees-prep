function particle = cellgap2(particle, swc, index_array,boundSize,lookup_table,step,perm_prob,crand)
    position = particle(1:3);
    state = particle(4);
    flag = particle(5);

    if (flag)
        flag = false;
    else
        % get positions
%         current = position;

        next = setnext(position', step,crand);

        % check positions
        inside = checkpos(next, swc, index_array, state,boundSize,lookup_table);

        if inside
            position = next;
            state = true;
        else
            % random walker was outside or crossed boundary
            % TODO ADD PROB back in... needs to be passed
            if (rand < perm_prob)
                position = next;
                state = false;
            else
                flag = true;
            end

        end

    end

    % rewrite particle
    particle(1:3) = position;
    particle(4) = state;
    particle(5) = flag;

    function next = setnext(position, step, crand)
%         v = normr(rand(1,3));
%         v = randn(1,1,3)';
%         v = v./sqrt(sum(v,2)); 
%         vect = random_unit_vector(3, 1);
%         vect = cos(180*pi*normr(rand(1,3)))';
        delta = crand * step;
        next = position + delta;

    end

    function inside = checkpos(next, swc, A, currstate,boundSize,lookup_table)
        % extract {child,parent} Ids of segments near location of randomwalker
        indicies = preprocesses(next,boundSize,lookup_table);

        % if there are elements of the cell near the random walker
        if (indicies ~= 0)

            % check if the random walker is inside these connections
            [currinside, nextinside] = check_connections(indicies, A, swc, next, currstate);

            % determine whether boundary was crossed
            % using states of current position and next
            inside = insideLogic(currinside, nextinside);
        else
            % if no elements are near the randomwalker then no boundaries were crossed
            inside = 0;
        end

    end

    function [currinside, nextinside] = check_connections(indicies, A, swc, next, currstate)
        % initialize current and next states
        currinside = currstate; nextinside = false;

        % extract child and parent ids
        children = A{indicies, 1}; parents = A{indicies, 2};

        % set xyz for inside-outside calculation
        nx0 = next(1); ny0 = next(2); nz0 = next(3);

        swc1 = swc(children, 2:5);
        swc2 = swc(parents, 2:5);

        % for each pair: check if point is inside
        for i = 1:length(children)

            % get base and target ids
            p1 = swc1(i, :);
            p2 = swc2(i, :);

            x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
            x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);

%             dist = memoized_distance(children(i));
            dist = (x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2;

            if r1 > r2
                tcn = pointbasedswc2v(nx0, ny0, nz0, x1, x2, y1, y2, z1, z2, r1, r2, false, dist);
            else
                tcn = pointbasedswc2v(nx0, ny0, nz0, x2, x1, y2, y1, z2, z1, r2, r1, false, dist);
            end

            % inside ith connection or already inside
            nextinside = tcn | nextinside;

            % if both positions are inside: break
            if currinside && nextinside
                break;
            end

        end

    end

    function inside = insideLogic(currinside, nextinside)
        % both positions within a connection
        if currinside && nextinside
            inside = 1;
            % disp("Was just inside")

            % current pos in, next pos out
        elseif currinside && ~nextinside
            % disp("next_outside")
            inside = 0;

            % current pos out, next pos in
        elseif ~currinside && nextinside

            disp("RW Outside NextInside: ");
            inside = 1;

            % both positions outside
        elseif ~currinside && ~nextinside
            disp("RW Outside: ");
            inside = 1;
        else
            disp("OUTSIDE: ")
            inside = 1;
        end

    end

    function indicies = preprocesses(next,boundSize,lookup_table)
            % convert float to voxel
            nvoxes = float2vox(next)';

            % logical array for elements within range 0 < n < boundSize
            nv = all(nvoxes > 0, 2) & all(nvoxes < boundSize, 2);

            % extract elements within range
            NVOXES = nvoxes(nv, :)';

            % convert subscript index to linear index
            % ^ (faster than using sub index to get values from lookup table)
            i1 = NVOXES(1,:);
            i2 = NVOXES(2,:);
            i3 = NVOXES(3,:);
            NVOXES = i1 + (i2-1)*boundSize(1) + (i3-1)*boundSize(1)*boundSize(2);
%             NVOXES = sub2ind(boundSize, NVOXES(1, :), NVOXES(2, :), NVOXES(3, :));

            % extract [child,parent] node ids
            nindicies = lookup_table(NVOXES);

            % return non-zero indicies or return 0
            indicies = combineind(nindicies');
    end

    function inds = combineind(nind)

        % extract non-zero elements
%         nidx = nind > 0; 
        tn = nind(nind>0);

        % logicals for non empty index arrays
        n = ~isempty(tn);

        if n
            % if there are non-zero elements return them
            inds = tn;
        else
            % return 0
            inds = 0;
        end

    end

end
