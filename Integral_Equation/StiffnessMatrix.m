function matA = StiffnessMatrix(vecMesh,a)
%Use to construct matrix A in the Galerkin method
intSize = length(vecMesh);
% vecDiag = zeros(intSize+1,1);
% vecSub = zeros(intSize+1,1);
% vecSup = vecSub;
matA = zeros(intSize,intSize);
intNumOfQuad = 2;

%fill thing in the middle
for i=2:intSize-2
    dblH = vecMesh(i) - vecMesh(i-1);
    dblH1 = vecMesh(i+1) - vecMesh(i);
    dblH2 = vecMesh(i+2) - vecMesh(i+1);
    
    %fill in supper diagonal
    %Find the integral of the first term
    [x1,w1] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
    funInt1 =@(x) ((-1/dblH1)*(x-vecMesh(i)) + 1).*(1/dblH1*(x-vecMesh(i)));
    dblInt1 = sum(feval(funInt1,x1).*w1);
    
    %Find the integral of the second terms - first term of K
    [y1,z1] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
    funInt2K1xy =@(x,y) a/pi*(1./((x-y).^2 + a^2)).*(1/dblH1*(y-vecMesh(i)));
    
    funInt2K1x =@(x) sum(feval(funInt2K1xy,x,y1).*z1);
    
    [y2,z2] = lgwt(intNumOfQuad,vecMesh(i+1),vecMesh(i+2));
    funInt2K2xy =@(x,y) a/pi*(1./((x-y).^2 + a^2)).*(-1/dblH2*(y-vecMesh(i+1)) + 1);
    
    funInt2K2x =@(x) sum(feval(funInt2K2xy,x,y2).*z2);
    
    
    [x2,w2] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
    
    funInt2 =@(x) 1/dblH*(x-vecMesh(i-1)).*(funInt2K1x(x) + funInt2K2x(x));
    dblInt2 = sum(feval(funInt2,x2).*w2);
    
    %Find the integral of the third terms - second term of K
    [x3,w3] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
    funInt3 =@(x) (-1/dblH1*(x-vecMesh(i)) + 1).*(funInt2K1x(x) + funInt2K2x(x));
    dblInt3 = sum(feval(funInt3,x3).*w3);
    
    matA(i,i+1) = dblInt1 - (dblInt2+dblInt3);
    
    %fill in the diagonal
    %Find integral of the first term
    [x1D1,w1D1] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
    funInt1D1 =@(x) (1/dblH*(x-vecMesh(i-1))).^2;
    [x1D2,w1D2] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
    funInt1D2 =@(x) (-1/dblH1*(x-vecMesh(i)) + 1).^2;
    dblInt1D = sum(feval(funInt1D1,x1D1).*w1D1) + sum(feval(funInt1D2,x1D2).*w1D2);
    
    %Find integral of the second term - first term of K
    [y1D,z1D] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
    funInt2K1Dxy =@(x,y) a/pi*(1./((x-y).^2 + a^2)).*(1/dblH*(y-vecMesh(i-1)));
    
    funInt2K1Dx =@(x) sum(feval(funInt2K1Dxy,x,y1D).*z1D);
    
    [y2D,z2D] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
    funInt2K2Dxy =@(x,y) a/pi*(1./((x-y).^2 + a^2)).*(-1/dblH1*(y-vecMesh(i)) + 1);
    
    funInt2K2Dx =@(x) sum(feval(funInt2K2Dxy,x,y2D).*z2D);
    
    [x2D,w2D] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
    
    funInt2D =@(x) 1/dblH*(x-vecMesh(i-1)).*(funInt2K1Dx(x) + funInt2K2Dx(x));
    dblInt2D = sum(feval(funInt2D,x2D).*w2D);
    
    %Find integral of the third term - second term of K
    [x3D,w3D] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
    funInt3D =@(x) (-1/dblH1*(x-vecMesh(i)) + 1).*(funInt2K1Dx(x) + funInt2K2Dx(x));
    dblInt3D = sum(feval(funInt3D,x3D).*w3D);
    
    matA(i,i) = dblInt1D - (dblInt2D+dblInt3D);
    
    %Fill in the sub diagonal
    %matA(i+1,i) = matA(i,i+1);
    %Find integral of the first term
    [x1L,w1L] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
    funInt1L =@(x) ((-1/dblH1)*(x-vecMesh(i)) + 1).*(1/dblH1*(x-vecMesh(i)));
    dblInt1L = sum(feval(funInt1L,x1L).*w1L);
    
    %Find integral of the second term - first term of K
    [y1L,z1L] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
    funInt2K1Lxy =@(x,y) a/pi*(1./((x-y).^2 + a^2)).*(1/dblH*(y-vecMesh(i-1)));
    
    funInt2K1Lx =@(x) sum(feval(funInt2K1Lxy,x,y1L).*z1L);
    
    [y2L,z2L] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
    funInt2K2Lxy =@(x,y) a/pi*(1./((x-y).^2 + a^2)).*(-1/dblH1*(y-vecMesh(i)) + 1);
    
    funInt2K2Lx =@(x) sum(feval(funInt2K2Lxy,x,y2L).*z2L);
    
    [x2L,w2L] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
    
    funInt2L =@(x) 1/dblH1*(x-vecMesh(i)).*(funInt2K1Lx(x) + funInt2K2Lx(x));
    dblInt2L = sum(feval(funInt2L,x2L).*w2L);
    
    %Find integral of the third term - second term of K
    [x3L,w3L] = lgwt(intNumOfQuad,vecMesh(i+1),vecMesh(i+2));
    funInt3L =@(x) (-1/dblH2*(x-vecMesh(i+1)) + 1).*(funInt2K1Lx(x) + funInt2K2Lx(x));
    dblInt3L = sum(feval(funInt3L,x3L).*w3L);
    
    matA(i+1,i) = dblInt1L - (dblInt2L+dblInt3L);
    
end
matA(intSize-1,intSize-1) = matA(intSize-2,intSize-2);

for i=2:intSize-1
   for j=2:intSize-1
      if matA(i,j) == 0
          
          dblHj = vecMesh(j) - vecMesh(j-1);
          dblHj1 = vecMesh(j+1) - vecMesh(j);
          dblHi = vecMesh(i) - vecMesh(i-1);
          dblHi1 = vecMesh(i+1) - vecMesh(i);
          
          %Find integral of the first term - first term of K
          [y1X,z1X] = lgwt(intNumOfQuad,vecMesh(j-1),vecMesh(j));
          funInt2K1Xxy =@(x,y) a/pi*(1./((x-y).^2 + a^2)).*(1/dblHj*(y-vecMesh(j-1)));

          funInt2K1Xx =@(x) sum(feval(funInt2K1Xxy,x,y1X).*z1X);

          [y2X,z2X] = lgwt(intNumOfQuad,vecMesh(j),vecMesh(j+1));
          funInt2K2Xxy =@(x,y) a/pi*(1./((x-y).^2 + a^2)).*(-1/dblHj1*(y-vecMesh(j)) + 1);

          funInt2K2Xx =@(x) sum(feval(funInt2K2Xxy,x,y2X).*z2X);

          [x2X,w2X] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));

          funInt2X =@(x) 1/dblHi*(x-vecMesh(i-1)).*(funInt2K1Xx(x) + funInt2K2Xx(x));
          dblInt2X = sum(feval(funInt2X,x2X).*w2X);

          %Find the integral of the second term - second term of K
          [x3X,w3X] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
          funInt3X =@(x) (-1/dblHi1*(x-vecMesh(i)) + 1).*(funInt2K1Xx(x) + funInt2K2Xx(x));
          dblInt3X = sum(feval(funInt3X,x3X).*w3X);
          
          matA(i,j) = -(dblInt2X + dblInt3X);
          
          
      end
   end
    
end

%Fill in the edge entries of the matrix

%First Row
i = 1;
for j=1:intSize
    dblInt1 = 0;
   if j == i+1
        %Find the integral of the first term
        dblH1 = vecMesh(i+1) - vecMesh(i);
        [x1,w1] = lgwt(intNumOfQuad,vecMesh(i),vecMesh(i+1));
        funInt1 =@(x) ((-1/dblH1)*(x-vecMesh(i)) + 1).*(1/dblH1*(x-vecMesh(i)));
        dblInt1 = sum(feval(funInt1,x1).*w1);   
   end
   dblInt2 = 0;
end


% %Fill the last entries of the matrix
% i = intSize-1;
% %fill in the diagonal
%     %Find integral of the first term
%     [x1D1,w1D1] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
%     funInt1D1 =@(x) (1/dblH*(x-vecMesh(i-1))).^2;
%     dblInt1D = sum(feval(funInt1D1,x1D1).*w1D1);
%     
%     %Find integral of the second term - first term of K
%     [y1D,z1D] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
%     funInt2K1Dxy =@(x,y) a/pi*(1./((x-y).^2 + a^2)).*(1/dblH*(y-vecMesh(i-1)));
%     
%     funInt2K1Dx =@(x) sum(feval(funInt2K1Dxy,x,y1D).*z1D);
%     
%     [x2D,w2D] = lgwt(intNumOfQuad,vecMesh(i-1),vecMesh(i));
%     funInt2D =@(x) 1/dblH*(x-vecMesh(i-1)).*(funInt2K1Dx(x));
%     dblInt2D = sum(feval(funInt2D,x2D).*w2D);
%     
%     matA(i,i) = dblInt1D - (dblInt2D+dblInt3D);

end

