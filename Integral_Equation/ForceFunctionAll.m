function [vecForce] = ForceFunctionAll(vecMesh,a)
%Calculate the forcing function
%funForce =@(x) 1 - 1/pi*(atan((1-x)/a) + atan((1+x)/a));
funForce = @(x) 1 ;%+ 1/pi*(-atan((x+1)/a) + atan((x-1)/a));
intSize = length(vecMesh);
vecForce = zeros(intSize,1);
intNumOfQuad =32;
for i=1:intSize
    switch i
        case 1
            %right
            dblH1 = vecMesh(i+1) - vecMesh(i);
            [x2,w2] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
            funInt2 =@(x) funForce(x).*(-1/dblH1*(x-vecMesh(i)) + 1);
            vecForce(i,1) = sum(feval(funInt2,x2).*w2);
        case intSize
            %left
            dblH = vecMesh(i) - vecMesh(i-1);
            [x1,w1] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
            funInt1 =@(x) funForce(x).*(1/dblH*(x-vecMesh(i-1)));
            vecForce(i,1) = sum(feval(funInt1,x1).*w1);
        otherwise
            %both
            dblH = vecMesh(i) - vecMesh(i-1);
            dblH1 = vecMesh(i+1) - vecMesh(i);
            %calculate the first integral
            [x1,w1] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
            funInt1 =@(x) funForce(x).*(1/dblH*(x-vecMesh(i-1)));
            [x2,w2] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
            funInt2 =@(x) funForce(x).*(-1/dblH1*(x-vecMesh(i)) + 1);
            vecForce(i,1) = sum(feval(funInt1,x1).*w1) + sum(feval(funInt2,x2).*w2);
    end
end

end

