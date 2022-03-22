function [vecExact] = exactfunction(vecMesh)
%EXACTFUNCTION Summary of this function goes here
%   Detailed explanation goes here
funExact = @(x) 1.9192 - 0.311717*x.^2 + 0.015676*x.^4 + 0.019682*x.^6 - 0.000373*x.^8;
vecExact = feval(funExact,vecMesh);
end

