% gpu test

A=ones(20,20,120);

x=A;
y=A;
x0=A;
y0=A;
zp=A;
lambda=A;

NA=A;
M=A;
ns=A;
ni=A;
DELTAt=A;
zd=A;
zdstar=A;
np=A;
na=A;

rho=A;
theta=A;
Bstar=A;
THETAplate=A;
THETAZernike=A; 

gpux=gpuArray(A);
gpuy=gpuArray(A);
gpux0=gpuArray(A);
gpuy0=gpuArray(A);
gpuzp=gpuArray(A);
gpulambda=gpuArray(A);

gpuNA=gpuArray(A);
gpuM=gpuArray(A);
gpuns=gpuArray(A);
gpuni=gpuArray(A);
gpuDELTAt=gpuArray(A);
gpuzd=gpuArray(A);
gpuzdstar=gpuArray(A);
gpunp=gpuArray(A);
gpuna=gpuArray(A);

gpurho=gpuArray(A);
gputheta=gpuArray(A);
gpuBstar=gpuArray(A);
gpuTHETAplate=gpuArray(A);
gpuTHETAZernike=gpuArray(A);

tic
Dnogpu=testfun(x, y, x0, y0, zp, lambda,...
                            NA, M, ns, ni, DELTAt, zd, zdstar, np, na,...
                            rho, theta, Bstar, THETAplate, THETAZernike);
nogpu=toc;

sprintf('no gpu = %0.3g',nogpu)

tic 
Dgpu=arrayfun(@testfun,gpux, gpuy, gpux0, gpuy0, gpuzp, gpulambda,...
            gpuNA, gpuM, gpuns, gpuni, gpuDELTAt, gpuzd, gpuzdstar, gpunp, gpuna,...
            gpurho, gputheta, gpuBstar, gpuTHETAplate, gpuTHETAZernike);
gpu=toc;

sprintf('gpu = %0.3g',gpu)


function out = testfun(x, y, x0, y0, zp, lambda,...
                            NA, M, ns, ni, DELTAt, zd, zdstar, np, na,...
                            rho, theta, Bstar, THETAplate, THETAZernike)
                        
 out= Bstar.*exp(1).^((sqrt(complex(-1))*2).*lambda.^(-1).*pi.*((-1).*DELTAt.*ni.*(1+(...
-1).*NA.^2.*ni.^(-2).*rho.^2).^(1/2)+(-1).*np.*THETAplate+(-1).*...
 na.*THETAZernike+(-1/2).*rho.*zd.^(-1).*zdstar.^(-1).*((-1).*zd+...
zdstar)+(-1).*ns.*(1+(-1).*NA.^2.*ns.^(-2).*rho.^2).^(1/2).*zp+...
M.^(-1).*NA.*rho.*((x+(-1).*x0).^2+(y+(-1).*y0).^2).^(1/2).*cos(...
theta+(-1).*atan2(x+(-1).*x0,y+(-1).*y0)))).*lambda.^(-1).*M.^(-2)...
.*NA.^2.*rho;

end
