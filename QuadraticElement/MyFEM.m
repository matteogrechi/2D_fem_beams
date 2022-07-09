classdef MyFEM
    properties
        beams  = {}; % Collection of beams objects
        constr = {}; % Constraints between nodes
        forces = {};
    end

    properties(Dependent)
        K; % Total stiffness
        f; % Generalized force vector
        A; % Constraint matrix
        b; % Constraint vector
        u; % Generalized nodes displacement
        l; % Lagrange multipliers
    end

    methods
        function obj = MyFEM(beams,constr,forces)
            arguments
                beams
                constr
                forces
            end

            obj.beams = beams;
            obj.constr = constr;
            obj.forces = forces;
        end

        function K = get.K(obj)
            nBeams = numel(obj.beams);
            K = sparse(9 * nBeams, 9 * nBeams);
            for i = 1:nBeams
                idx = 9 * (i - 1) + 1;
                K(idx:idx+8,idx:idx+8) = obj.beams{i}.K;
            end
        end

        function f = get.f(obj)
            flen = numel(obj.beams) * 9;
            f = sparse(flen,1);
            % Forces contribution
            for i = 1:numel(obj.forces)
                f = f + obj.forces{i}.f;
            end
            % Beams weight contribution
            for i = 1:numel(obj.beams)
                idx = 9 * (i - 1) + 1;
                f(idx:idx+8) = f(idx:idx+8) + obj.beams{i}.w;
            end
        end

        function A = get.A(obj)
            Acols = size(obj.K,2);
            A = sparse(0,Acols);
            for i = 1:numel(obj.constr)
                A = [A; obj.constr{i}.A];
            end
        end

        function b = get.b(obj)
            b = sparse(0,1);
            for i = 1:numel(obj.constr)
                b = [b; obj.constr{i}.b];
            end
        end

        function u = get.u(obj)
            U = sparse(null(full(obj.A)));
            u_ = obj.A\obj.b;

            lambda = (U' * obj.K * U) \ (U' * obj.f - U' * obj.K * u_);
            u = u_ + U * lambda;
        end

        function l = get.l(obj)
            l = (obj.A * obj.A') \ (obj.A * obj.f - obj.A * obj.K * obj.u);
        end

        function [u,w,theta] = getNodeDisplacement(obj,beamIdx,nodeIdx)
            idx = 9 * (beamIdx - 1) + 3 * (nodeIdx - 1) + 1;
            
            u = full(obj.u(idx));
            w = full(obj.u(idx + 1));
            theta = full(obj.u(idx + 2));
        end
    end
end