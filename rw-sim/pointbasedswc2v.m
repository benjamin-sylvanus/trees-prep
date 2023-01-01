% Replaced Calulation for r comparison
function pos = pointbasedswc2v(x0, y0, z0, x1, x2, y1, y2, z1, z2, r1, r2, emptylogical)

    if r2 > r1
        pos = pointbasedswc2v(x0, y0, z0, x2, x1, y2, y1, z2, z1, r2, r1, emptylogical);
    else

        t = ((x0 - x1) * (x2 - x1) + (y0 - y1) * (y2 - y1) + (z0 - z1) * (z2 - z1)) ./ ...
            ((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);
        x = x1 + (x2 - x1) * t;
        y = y1 + (y2 - y1) * t;
        z = z1 + (z2 - z1) * t;

        if (x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2 < r1 ^ 2
            list1 = emptylogical;
            list2 = ~list1;
        else
            list1 = (x - x1) .* (x - x2) + (y - y1) .* (y - y2) + (z - z1) .* (z - z2) < 0;
            list2 = ~list1;

        end

        dist2 = (x0 - x) .^ 2 + (y0 - y) .^ 2 + (z0 - z) .^ 2;

        if list1
            rd = abs(r1 - r2);

            % distance from orthogonal vector to p2
            l = sqrt((x - x2) .^ 2 + (y - y2) .^ 2 + (z - z2) .^ 2);

            % distance from p1 -> p2
            L = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);

            c = (rd * l) ./ L;
            r = (c + r2) ./ sqrt(1 - ((rd / L) .^ 2));
            pos1 = dist2 < (r .^ 2); % smaller in one line and less than and equal
            pos = pos1;
        else
            pos2 = (((x0 - x1) .^ 2 + (y0 - y1) .^ 2 + (z0 - z1) .^ 2) < ((r1) ^ 2)) | ...
                (((x0 - x2) .^ 2 + (y0 - y2) .^ 2 + (z0 - z2) .^ 2) < ((r2) ^ 2));
            pos = pos2;
        end

        %     r = r1 + sqrt((x-x1).^2 + (y-y1).^2 + (z-z1).^2) /...
        %         sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) * (r2-r1);

        %     r = ( c + r2 ) / (sqrt ( 1 - ( |r1-r2 | / l ) )

        %     c = ( |r1 - r2| * l ) / L

        % rd = abs(r1 - r2);

        % % distance from orthogonal vector to p2
        % l = sqrt((x - x2) .^ 2 + (y - y2) .^ 2 + (z - z2) .^ 2);

        % % distance from p1 -> p2
        % L = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);

        % c = (rd * l) ./ L;
        % r = (c + r2) ./ sqrt(1 - ((rd / L) .^ 2));

        % pos1 = dist2 < (r .^ 2); % smaller in one line and less than and equal
        % pos2 = (((x0 - x1) .^ 2 + (y0 - y1) .^ 2 + (z0 - z1) .^ 2) < ((r1) ^ 2)) | ...
        %     (((x0 - x2) .^ 2 + (y0 - y2) .^ 2 + (z0 - z2) .^ 2) < ((r2) ^ 2));

        % % faster
        % pos = bitor(bitand(list1, pos1), bitand(list2, pos2));
        % % pos(list1) = pos1(list1);
        % % pos(list2) = pos2(list2);
    end

end
