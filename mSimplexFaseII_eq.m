function[xo,zo, ban, iter,B] = mSimplexFaseII_eq(A, b, c, jB)
% minimizar c^T x
% sujeto a Ax = b , x >= 0
%
% In : A ... mxn matrix
% b ... vector columna con m renglones
% c ... vector columna con n renglones
% jB ... vector de indices de la SBF inicial.
%
% Out: xo ... SFB ´optima del problema
% zo ... valor ´optimo del problema
% ban ... indica casos:
% 0 ... si se encontro una soluci´on ´optima
% 1 ... si la funci´on objectivo no es acotada.
% iter ... es el numer´o de iteraciones (cambios de variables basicas)
% que hizo el m´etodo
% jB ... vector de indices de la SBF optima

[m,n] = size(A);
iter = 0;
ban = 0;
stop = false;

%jB = sort(jB);
jN = setdiff(1:n,jB)';
tab = zeros(m+1,n+1);

%CONSTRUIR TABLEAU - proposición 2.1
Ab = A(:,jB);
An = A(:,jN);
cb = c(jB);
cn = c(jN);
H = Ab\An;
h = Ab\b;
alp = cb'*h;
rn = cb'*H - cn';
tab(1:end-1,jB) = eye(m);
tab(1:end-1,jN) = H;
tab(1:end-1,end) = h(:);
maxiter=1000;
tab(end, end) = alp;
tab(end,jN) = rn';
%disp(tab);

while (~stop && iter <maxiter) 
    if any(tab(end,1:n) > 0) 
        %disp(tab)
        [~,e] = max(tab(end,1:n)); %valor de xe, indice de xe
        ind_e = find(jN==e);
        
%     Usando R de max descenso
        if all(tab(1:m,e) <= 0 )
            ban = 1;
            %STOP, no acotada
            break
        else
            X=tab(1:m,end);
            Y=tab(1:m,e);            
            aux=X./Y;
            min_init = false;
            ent_posit = false;
            for k = 1:size(aux)
                if aux(k)>0
                    if ~min_init
                        minv=aux(k);
                        min_init = true;
                        s = k;
                        ent_posit = true;
                    else
                        if aux(k)< minv
                            minv = aux(k);
                            s = k;
                            ent_posit = true;
                        end
                    end
                end
            end        
            if ~ent_posit
                ban = 1;
                break;
            end
            
            
            %normalizar
            tab(s,1:end) = tab(s,1:end)./tab(s,e);                 
            for i = 1:m+1
                if i ~= s %hace a los otros renglones cero
                    tab(i,:) = tab(i,:) - tab(i,e)*tab(s,:);
                end
            end
            T = ["Iteración número: ", num2str(iter +1)];
            %disp(T)
            %disp(tab)
            iter = iter + 1;
            aux = jN(ind_e);
            
            %size(jB)
            %size(jN)
            jN(ind_e) = jB(s);
            jB(s) = aux;
            
        end  
    else
        break
    end     
end    

xo_aux = tab(1:m, n+1);
m_aux = max(m,n-m);
xo = zeros(m_aux,1);
B = jB;
for i = 1:m
    if B(i)<m_aux
        xo(B(i))=xo_aux(i);
    end
end

    
        

zo = tab(end, end);

end
