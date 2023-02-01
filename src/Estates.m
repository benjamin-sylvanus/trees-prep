classdef Estates < handle

    properties
        enum;
    end
    
    methods
        function obj = Estates(curr,next)
            if (curr&&next)
                obj.enum=enumstates.inside;
            end
             if (~curr&&~next)
                 obj.enum=enumstates.outside;
             end
            if (curr&&~next)
                obj.enum=enumstates.crossout;
            end
            if (~curr&&next)
                obj.enum =  enumstates.crossin;
            end
        end
    end

end



