function [matA] = StiffnessMatrixAll(vecMesh,a)
%Use to compute the matrix A
intSize = length(vecMesh);
matA = zeros(intSize,intSize);
intNumOfQuad = 32;


for i = 1:intSize
   for j= 1:intSize
       %Find the first part of the integral
       dblInt1 = 0;
       if i == j%diagonal
           dblInt1D = 0;
           if i == 1
               %right
               dblH1 = vecMesh(i+1) - vecMesh(i);
               [x1D2,w1D2] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
               funInt1D2 =@(x) (-1/dblH1*(x-vecMesh(i)) + 1).^2;
               dblInt1D = sum(feval(funInt1D2,x1D2).*w1D2);
           elseif i == intSize
               %left
               dblH = vecMesh(i) - vecMesh(i-1);
               [x1D1,w1D1] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
               funInt1D1 =@(x) (1/dblH*(x-vecMesh(i-1))).^2;
               dblInt1D = sum(feval(funInt1D1,x1D1).*w1D1);
                   
           else
               %both
               dblH = vecMesh(i) - vecMesh(i-1);
               dblH1 = vecMesh(i+1) - vecMesh(i);
               
               [x1D1,w1D1] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
               funInt1D1 =@(x) (1/dblH*(x-vecMesh(i-1))).^2;
               [x1D2,w1D2] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
               funInt1D2 =@(x) (-1/dblH1*(x-vecMesh(i)) + 1).^2;
               dblInt1D = sum(feval(funInt1D1,x1D1).*w1D1) + sum(feval(funInt1D2,x1D2).*w1D2);
           end
           dblInt1 = dblInt1D;
           
       elseif abs(i-j) == 1%the upper and lower diagonal, this value is the same
           k=min(i,j);
           %right of min, left of max
           dblH1 = vecMesh(k+1) - vecMesh(k);
           [x1,w1] = lgwt(intNumOfQuad,vecMesh(k),vecMesh(k+1));
           funInt1 =@(x) ((-1/dblH1)*(x-vecMesh(k)) + 1).*(1/dblH1*(x-vecMesh(k)));
           dblInt1 = sum(feval(funInt1,x1).*w1);
           %test = integral(funInt1,vecMesh(k),vecMesh(k+1));

       end
       
       %Find the second integral
       %Find the function for the outer part
       %Dummy intialization
       funInt2OutP1 =@(x) 0;
       funInt2OutP2 =@(x) 0;
       [x2P1,w2P1] =lgwt(intNumOfQuad,0,1);
       x2P2 = x2P1;
       w2P2 = w2P1;
       switch i
           case 1
           %right
           [x2P2,w2P2] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
           funInt2OutP2 =@(x) -1/(vecMesh(i+1)-vecMesh(i))*(x-vecMesh(i)) + 1;
           case intSize
           %left
           [x2P1,w2P1] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
           funInt2OutP1=@(x) 1/(vecMesh(i)-vecMesh(i-1))*(x-vecMesh(i-1));
           otherwise
           %both
           [x2P1,w2P1] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
           funInt2OutP1=@(x) 1/(vecMesh(i)-vecMesh(i-1))*(x-vecMesh(i-1));
           
           [x2P2,w2P2] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
           funInt2OutP2 =@(x) -1/(vecMesh(i+1)-vecMesh(i))*(x-vecMesh(i)) + 1;
       end
       
       %Find the function for the inner part
       %Dummy initialization
       funInt2In =@(x) 0;
       
       [y1P1,z1P1] =lgwt(intNumOfQuad,0,1);
       y1P2 = y1P1;
       z1P2 = z1P1;
       
       switch j
           case 1
               %right
               [y1P2,z1P2] = lgwt(intNumOfQuad,vecMesh(j),vecMesh(j+1));
               funInt2InP2xy =@(x,y) a/pi*(1./((x-y).^2 + a^2))...
                   .*(-1/(vecMesh(j+1)-vecMesh(j)).*(y-vecMesh(j)) + 1);
               funInt2InP2 =@(x) sum(feval(funInt2InP2xy,x,y1P2).*z1P2);
               funInt2In =@(x) funInt2InP2(x);
           case intSize
               %left
               [y1P1,z1P1] = lgwt(intNumOfQuad,vecMesh(j-1),vecMesh(j));
               funInt2InP1xy =@(x,y) a/pi*(1./((x-y).^2 + a^2))...
                   .*(1/(vecMesh(j)-vecMesh(j-1)).*(y-vecMesh(j-1)));
               funInt2InP1 =@(x) sum(feval(funInt2InP1xy,x,y1P1).*z1P1);
               funInt2In =@(x) funInt2InP1(x);
               
           otherwise
               %left
               [y1P1,z1P1] = lgwt(intNumOfQuad,vecMesh(j-1),vecMesh(j));
               funInt2InP1xy =@(x,y) a/pi*(1./((x-y).^2 + a^2))...
                   .*(1/(vecMesh(j)-vecMesh(j-1)).*(y-vecMesh(j-1)));
               funInt2InP1 =@(x) sum(feval(funInt2InP1xy,x,y1P1).*z1P1);
               
               %right
               [y1P2,z1P2] = lgwt(intNumOfQuad,vecMesh(j),vecMesh(j+1));
               funInt2InP2xy =@(x,y) a/pi*(1./((x-y).^2 + a^2))...
                   .*(-1/(vecMesh(j+1)-vecMesh(j)).*(y-vecMesh(j)) + 1);
               funInt2InP2 =@(x) sum(feval(funInt2InP2xy,x,y1P2).*z1P2);
               
               funInt2In =@(x) funInt2InP1(x) + funInt2InP2(x);
               
       end
           
       %value of the second integral
       funInt2P1 =@(x) funInt2OutP1(x).*funInt2In(x);
       funInt2P2 =@(x) funInt2OutP2(x).*funInt2In(x);
       
       dblInt2 = sum(feval(funInt2P1,x2P1).*w2P1) + sum(feval(funInt2P2,x2P2).*w2P2);
       
       if i ==2 && j==3
          test = 1; 
       end
       
       matA(i,j) = dblInt1 - dblInt2;
       
   end    
end

end

