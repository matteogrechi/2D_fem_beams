function K = getKloc(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
E = obj.E;
A = obj.A;
J = obj.J;
L = obj.l;

Ku = E*A/L * [  7/3, -8/3,  1/3 ;
               -8/3, 16/3, -8/3 ;
                1/3, -8/3,  7/3 ];

Kw = E*J/L^3 * [     5092/35,  (1138*L)/35,     -512/5,   (384*L)/7,     -1508/35,   (242*L)/35 ;
                 (1138*L)/35, (332*L^2)/35, -(128*L)/5,  (64*L^2)/7,  -(242*L)/35,  (38*L^2)/35 ;
                      -512/5,   -(128*L)/5,     1024/5,           0,       -512/5,    (128*L)/5 ;
                   (384*L)/7,   (64*L^2)/7,          0, (256*L^2)/7,   -(384*L)/7,   (64*L^2)/7 ;
                    -1508/35,  -(242*L)/35,     -512/5,  -(384*L)/7,      5092/35, -(1138*L)/35 ;
                  (242*L)/35,  (38*L^2)/35,  (128*L)/5,  (64*L^2)/7, -(1138*L)/35, (332*L^2)/35 ];

K = zeros(9);
for i = 1:3
    for j = 1:3
        rw = 2 * (i-1) + 1;
        cw = 2 * (j-1) + 1;
        r = 3 * (i-1) + 1;
        c = 3 * (j-1) + 1;
        
        K(r,c) = Ku(i,j);
        K(r+1:r+2,c+1:c+2) = Kw(rw:rw+1,cw:cw+1);
    end
end

end