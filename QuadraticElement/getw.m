function w = getw(obj)
%GETW Summary of this function goes here
%   Detailed explanation goes here
g = 9.8086;
theta = obj.theta;
L = obj.l;
m = obj.rho * (obj.A * 1e-6) * (L * 1e-3);

w = [          (m*g*sin(2*theta))/30 ;
      -(m*g*(2*cos(theta)^2 + 5))/30 ;
              -(L*m*g*cos(theta))/60 ;
              -(m*g*sin(2*theta))/15 ;
       (2*m*g*(cos(theta)^2 - 5))/15 ;
                                   0 ;
               (m*g*sin(2*theta))/30 ;
      -(m*g*(2*cos(theta)^2 + 5))/30 ;
               (L*m*g*cos(theta))/60 ];

end

