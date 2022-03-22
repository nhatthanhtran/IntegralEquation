function funKernel = kernel(a)
%Love kernel
funKernel=@(x,y) a/pi*(1/((x-y).^2) + a^2);
end

