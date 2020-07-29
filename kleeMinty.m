%Script empírico KleeMinty
t = zeros(1,8);
iters = zeros(1,8);
tama= zeros(1,8);
for i = 1:8
    tic;
    m = i+2
    [c, A, b] = generaKleeMinty(m);
    [xo, zo, ban, iter,B] = mSimplex_leq(A, b, c);
    tama(i)=m
    iters(i) = iter
    t(i) = toc
end

figure(1);
plot(tama,iters,'r','LineWidth',2);
title('# iteraciones')

figure(2);
plot(tama,t,'b','LineWidth',2);
title('tiempo')