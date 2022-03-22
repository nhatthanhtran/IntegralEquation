clc;clear all;
close all;
%Initialize

intNumOfTrials = 5;
a=1;

matError = zeros(intNumOfTrials,3);
for i=1:intNumOfTrials
    
    intGridpoints = 2^i*10; 
    dx = 2/intGridpoints;
    vecMesh = -1:dx:1;
    vecExact = ones(length(vecMesh),1);
    matA = StiffnessMatrixAll(vecMesh,a);
    vecForce = ForceFunctionAll(vecMesh,a);

    vecU = matA\vecForce;
    matError(i,1) = dx;
    matError(i,2) = norm(vecU - vecExact,2)^2;
    %matError(i,2) = sum((vecU - vecExact).^2);
    if i ~= 1
        matError(i,3) = log(matError(i,2)/matError(i-1,2))/log(matError(i,1)/matError(i-1,1));
    end
    
end

%Calculate ROC

disp(matError)



