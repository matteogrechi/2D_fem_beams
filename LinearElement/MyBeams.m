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
            RR = blkdiag(obj.R,obj.R);

            w = RR * obj.w;
        end

        function w = get.w(obj)
            mg = obj.rho * (obj.A * 1e-6) * (obj.l * 1e-3) * 9.8086;
            
            w = zeros(6,1);
            w(1) = 0;
            w(2) = -mg/2;
            w(3) = -mg*obj.l*cos(obj.theta)/12;
            w(4) = 0;
            w(5) = -mg/2;
            w(6) = mg*obj.l*cos(obj.theta)/12;
        end

        function K = get.Kloc(obj)
            E_ = obj.E;
            l_ = obj.l;
            A_ = obj.A;
            J_ = obj.J;
            
            Ku = E_*A_/l_ * [  1 -1 ; 
                              -1  1 ];

            Kw_11 = E_*J_/l_^3 * [   12    6*l_ ;
                                   6*l_ 4*l_*l_ ];

            Kw_12 = E_*J_/l_^3 * [   -12    6*l_ ;
                                   -6*l_ 2*l_*l_ ];
            
            Kw_22 = E_*J_/l_^3 * [    12   -6*l_ ;
                                   -6*l_ 4*l_*l_ ];

            K = zeros(6,6);

            K(1:3,1:3) = [    Ku(1,1) zeros(1,2) ;
                           zeros(2,1)      Kw_11 ];

            K(1:3,4:6) = [    Ku(1,2) zeros(1,2) ;
                           zeros(2,1)      Kw_12 ];

            K(4:6,4:6) = [    Ku(2,2) zeros(1,2) ;
                           zeros(2,1)      Kw_22 ];

            K(4:6,1:3) = K(1:3,4:6)';
        end

        function K = get.K(obj)
            RR = blkdiag(obj.R,obj.R);

            K =  RR'*obj.Kloc*RR;
        end
    end
end

