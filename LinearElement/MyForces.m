classdef MyForces
    %MYFORCES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nBeams;
        beamIdx;
        nodeIdx;
        Fx; % [N]
        Fy; % [N]
        Mz; % [N m]
    end

    properties(Dependent)
        f; % Generalized force vector
    end
    
    methods
        function obj = MyForces(nBeams,beamIdx,nodeIdx,Fx,Fy,Mz)
            arguments
                nBeams {mustBeInteger,mustBePositive}
                beamIdx {mustBeInteger,mustBePositive,mustBeLessThanOrEqual(beamIdx,nBeams)}
                nodeIdx {mustBeMember(nodeIdx,[1,2])}
                Fx double = 0
                Fy double = 0
                Mz double = 0
            end
            obj.nBeams = nBeams;
            obj.beamIdx = beamIdx;
            obj.nodeIdx = nodeIdx;
            obj.Fx = Fx;
            obj.Fy = Fy;
            obj.Mz = Mz;
        end
        
        function f = get.f(obj)
            flen = obj.nBeams * 2 * 3;
            idx = 6 * (obj.beamIdx - 1) + 3 * (obj.nodeIdx - 1) + 1;
            
            f = sparse(flen,1);
            f(    idx) = obj.Fx;
            f(idx + 1) = obj.Fy;
            f(idx + 2) = obj.Mz * 1e3; % From N*m to N*mm
        end
    end
end

