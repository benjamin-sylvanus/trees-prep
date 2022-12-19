classdef randomwalker
    %RANDOMWALKER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        curr;
        next;
        step;
    end
    
    methods
        function obj = randomwalker(currentPosition,step)
            %RANDOMWALKER Construct an instance of this class
            %   Detailed explanation goes here
            obj.curr = currentPosition;
        end

        function obj = initPosition(obj,flag,sx,sy,sz)

            if flag
%                 init inside cell 
%                 while loop until a point is inside the volume
                  
            else

                x = randi(sx,1);
                y = randi(sy,1);
                z = randi(sz,1);
                obj.curr = [x;y;z];

            end
        end

        function obj = setnext(obj)
            vect = random_unit_vector(3,1);
            delta = vect * obj.step;
            pos2 = obj.curr + delta;
            obj.next = pos2;
        end
        
    end
end

