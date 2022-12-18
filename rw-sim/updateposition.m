function [pos, pos2] = updateposition(step, pos)
    vect = random_unit_vector(3, 1);
    delta = vect * step;
    pos2 = pos + delta;
end
