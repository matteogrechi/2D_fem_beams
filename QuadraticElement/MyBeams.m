classdef MyBeams
    %MYBEAMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        E;     % Young modulus [MPa]
        rho;   % Material density [kg/m^3]
        A;     % Section area of the beam [mm^2]
        J;     % Moment of inertia of the beam [mm^4]
        l;     % Length of the beam [mm]
        theta; % Beam angle [rad]
    end

    properties(Dependent)
        R;    % Rotation matrix from global to local frame
        wloc; % Weight force in local frame
        w;    % Weight force in global frame (g=9.8086m/s^2 in -y global)
        Kloc; % Local frame stiffness matrix
        K;    % Global frame stiffness matrix
    end
    
    methods
        function obj = MyBeams(E,rho,A,J,l,theta)
            arguments
                E
                rho
                A
                J
                l
                theta = 0;
            end

            obj.E = E;
            obj.rho = rho;
            obj.A = A;
            obj.J = J;
            obj.l = l;
            obj.theta = theta;
        end

        function R = get.R(obj)
            th = obj.theta;

            R = [ cos(th) sin(th) 0;
                 -sin(th) cos(th) 0;
                        0       0 1];
        end

        function w = get.wloc(obj)
            RR = blkdiag(obj.R,obj.R,obj.R);

            w = RR * obj.w;
        end

        function w = get.w(obj)
            w = getw(obj);
        end

        function K = get.Kloc(obj)
            K = getKloc(obj);
        end

        function K = get.K(obj)
            RR = blkdiag(obj.R,obj.R,obj.R);

            K =  RR'*obj.Kloc*RR;
            
            % For numerical reasons the matrix is non symmetric
            K = 0.5 * K + 0.5 * K';
        end
    end
end

