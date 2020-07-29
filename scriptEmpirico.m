% % generar dimensiones del problema
 m = round(10*exp(log(20)*rand()));
 n = round(10*exp(log(20)*rand()));
% 
% % generar A, b, c
 sigma = 100;
 A = round(sigma*randn(m,n));
 b = round(sigma*abs(randn(m,1)));
 c = round(sigma*randn(n,1));

res = []; %Matriz de resultados
for i=1:50
% % generar dimensiones del problema
 m = round(10*exp(log(20)*rand()));
 n = round(10*exp(log(20)*rand()));
% 
% % generar A, b, c
 sigma = 100;
 A = round(sigma*randn(m,n));
 b = round(sigma*abs(randn(m,1)));
 c = round(sigma*randn(n,1));
 [xo, zo, ban, iter,B] = mSimplex_leq(A, b, c);
 res=[res;iter,m,n,ban];
 
end


bounded = [];
unbounded = [];
%A matrtiz 50x4 - col(1:3)={#iter,m,n} col(4) = {acotado, no acotado}
for i=1:size(res(:,1))
    if res(i,4)==1
        bounded = [bounded;[res(i,1),min(res(i,[2,3]))]];
    else
        unbounded = [unbounded;[res(i,1),min(res(i,[2,3]))]];
    end
end
bounded;
unbounded;


scatter(bounded(:,2),bounded(:,1),'b', 'filled','Marker','s')
hold on
scatter(unbounded(:,2), unbounded(:,1),'r','filled')
hold off
xlabel('min(m,n)', 'fontsize', 14);
ylabel('#it','fontsize', 14);
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca, 'XMinorTick','on')
set(gca, 'YMinorTick', 'on')
legend('Bounded','Unbounded')
grid on