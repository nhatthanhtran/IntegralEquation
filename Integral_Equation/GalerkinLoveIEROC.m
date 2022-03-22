clc;clear all;
close all;
format short e;
%Initialize
dx = 0.01;
vecMesh = -1:dx:1;
vecExact = ones(length(vecMesh),1);
% intNumData =50;
% Domain= [-1 1];
% vecMesh = (Domain(1,2) - Domain(1,1))*rand(1,intNumData) + Domain(1,1);
% vecMesh = sort(vecMesh);
% a = 0.5;
% 
% matA = StiffnessMatrixAll(vecMesh,a);
% vecForce = ForceFunctionAll(vecMesh,a);
% 
% %matA =matA(2:end-1,2:end-1);
% %vecForce = vecForce(2:end-1,1);
% 
% vecU = matA\vecForce;
% % 
% % plot(vecMesh,vecExact,vecMesh,vecU);
% % legend('exact', 'approx')
% % axis([-1 1 0.9 1.1])
% % xlabel('x')
% % ylabel('u(x)')
% vecExactTest = exactfunction(vecMesh);
% vecExactTestCheb = exactfunctionCheb(vecMesh);
% plot(vecMesh,vecU,vecMesh,vecExactTestCheb',vecMesh,vecExactTest)
% legend('approx Galerkin','approx Chebyshev', 'approx iterative')
% %axis([-1 1 0 2.1])
% xlabel('x')
% ylabel('u(x)')

% figure(2)
% plot(vecMesh,vecU)
% % funForce = @(x) 1 - 1/pi*(atan((x+1)/a) - atan((x-1)/a));
% % test = feval(funForce,vecMesh);
% AbsError = abs(vecU-vecExactTestCheb')';
vecA = [1 0.5 0.25 0.125 0.0125];
figure(3)
for i = 1:length(vecA)
    matA = StiffnessMatrixAll(vecMesh,vecA(i));
    vecForce = ForceFunctionAll(vecMesh,vecA(i));
    vecU = matA\vecForce;
    plot(vecMesh,vecU)
    hold on
end
legend('1','0.5','0.25','0.125','0.0125')
xlabel('x')
ylabel('u(x)')

