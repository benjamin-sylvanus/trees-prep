classdef enumstates

    properties
        cstate;
        nstate;
        state; 
    end

    methods
        function obj = enumstates(c,n,str)
            obj.cstate = c;
            obj.nstate = n;
            obj.state = str;
        end
    end
        
    enumeration
        inside (1,1, "INSIDE");
        outside (0,0,"OUTSIDE");
        crossin (0,1, "CROSSIN");
        crossout (1,0, "CROSSOUT");
    end
end

