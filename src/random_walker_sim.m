classdef random_walker_sim
    %RANDOM_WALKER_SIM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        lookup_table;
        index_array;
        pairs;
        swc;
        step_size;
        randomwalker;
    end

    methods (Access = public)

        function obj = random_walker_sim(inputArg1, inputArg2)
            %RANDOM_WALKER_SIM Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

    end

end
