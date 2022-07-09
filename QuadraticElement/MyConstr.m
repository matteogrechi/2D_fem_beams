classdef MyConstr
    %MYCONSTR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nBeams;   % Number of beams
        beam1idx; % Beam 1 idx
        node1idx; % 1 or 2 depending on which node on the beam is considered
        beam2idx; % Beam 2 idx
        node2idx; % same as node1idx
        type;     % constraint type
    end

    properties(Dependent)
        A; % A matrix (used in augmented system with lagrange multipliers)
        b; % b vector (used in augmented system with lagrange multipliers)
    end
    
    methods
        function obj = MyConstr(nBeams,type,beam1idx,node1idx,beam2idx,node2idx)
            arguments
                nBeams {mustBeInteger,mustBePositive}
                type {mustBeMember(type,["clamp","hinge"])}
                beam1idx {mustBeInteger,mustBePositive,mustBeLessThanOrEqual(beam1idx,nBeams)}
                node1idx {mustBeMember(node1idx,[1,2,3])}
                beam2idx {mustBeInteger,mustBeNonnegative,mustBeLessThanOrEqual(beam2idx,nBeams)} = 0
                node2idx {mustBeMember(node2idx,[1,2,3])} = 1
            end

            obj.nBeams = nBeams;
            obj.type = type;
            obj.beam1idx = beam1idx;
            obj.node1idx = node1idx;
            obj.beam2idx = beam2idx;
            obj.node2idx = node2idx;
        end
        
        function A = get.A(obj)
            Acols = 9 * obj.nBeams;
            idx1 = 9 * (obj.beam1idx - 1) + 3 * (obj.node1idx - 1) + 1;
            idx2 = 9 * (obj.beam2idx - 1) + 3 * (obj.node2idx - 1) + 1;
            % idx2 is < 0 when the contraint is to the ground
            idx2 = max(idx2,0);

            switch obj.type
                case "hinge"
                    A = sparse(2,Acols);
                    A(1,  idx1) = 1;
                    A(2,idx1+1) = 1;
                    if(idx2 > 0)
                        A(1,  idx2) = -1;
                        A(2,idx2+1) = -1;
                    end
               
                case "clamp"
                    A = sparse(3,Acols);
                    A(1,  idx1) = 1;
                    A(2,idx1+1) = 1;
                    A(3,idx1+2) = 1;
                    if(idx2 > 0)
                        A(1,  idx2) = -1;
                        A(2,idx2+1) = -1;
                        A(3,idx2+2) = -1;
                    end
            end 
        end
    
        function b = get.b(obj)
            switch obj.type
                case "hinge"
                    b = sparse(2,1);

                case "clamp"
                    b = sparse(3,1);
            end
        end
    end
end

