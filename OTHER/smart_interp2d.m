function out = smart_interp2d(A, ydim, zdim)


ydim0 = size(A,1);
zdim0 = size(A,2);

ylin = linspace(-1,1,ydim0);
zlin = linspace(-1,1,zdim0);

ylinq = linspace(-1,1,ydim);
zlinq = linspace(-1,1,zdim);



[Z, Y] = meshgrid(zlin, ylin);
[Zq, Yq] = meshgrid(zlinq, ylinq);



out = interp2(Z, Y, A+10, Zq, Yq)-10;
