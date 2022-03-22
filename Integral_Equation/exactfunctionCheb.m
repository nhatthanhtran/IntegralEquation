function [vecExact] = exactfunctionCheb(vecMesh)
%EXACTFUNCTIONCHEV Summary of this function goes here
%   Detailed explanation goes here
funExact =@(x) 3.5488889/2 -0.1400431*cos(2*acos(x)) + 0.0049620*cos(4*acos(x))...
    +0.0003763*cos(6*acos(x)) - 0.0000432*cos(8*acos(x)) - 0.0000016*cos(10*acos(x));
vecExact = feval(funExact,vecMesh);
end

