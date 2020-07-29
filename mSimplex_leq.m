function[xo, zo, ban, iter,B] =mSimplex_leq(A, b, c)
% purpose: Versión del Simplex (Fase I, wrap)
% minimizar c^T x
% sujeto a Ax <= b , x >= 0
%
% In : A ... mxn matrix
% b ... vector columna con m renglones
% c ... vector columna con n renglones
%
% Out: xo ... SFB óptima del problema
% zo ... valor óptimo del problema
% ban ... indica casos:
% -1 ... si el conjunto factible es vacio
% 0 ... si se encontro una solución óptima
% 1 ... si la función objectivo no es acotada.
% iter ... es el numeró de iteraciones (cambios de variables basicas)
% que hizo el m´ etodo
% B ... vector de indices de la SBF optima
%

[m,n] = size(A);
iter = 1;
ban = 0;
c_orig = [c;zeros(m,1)];

if(b>=0) %case1: 0 pertenece a Cf
    A = [A,eye(m)];
    jB=(n+1:m+n)';
    c = [c;zeros(m,1)];
    [xo, zo, ban, iter, B] = mSimplexFaseII_eq(A,b,c,jB);
    xo = xo(1:n,1);
else %Hay entradas b(i) < 0
    %Contruir tableau del prob auxiliar
    jB = (n+1:n+m)';
    jN = [1:n,n+m+1]';
    
    tab = zeros(m+1,n+m+2);
    tab(1:m,1:n) = A;
    tab(1:end,end-1) = -ones(m+1,1);
    tab(1:m,end) = b(:);
    tab(1:m,n+1:n+m) = eye(m);
    
    %entra x_art sale xs
    [~,ind_s]=min(tab(:,n+m+2));
    var_e = n+m+1;
    ind_e = find(jN==var_e);
    
    aux = jN(ind_e);
    jN(ind_e) = jB(ind_s);
    jB(ind_s) = aux;
    
    tab(ind_s,1:end) = -tab(ind_s,1:end);
    
    for i = 1:m+1
        if i ~= ind_s %hace a los otros renglones cero
            tab(i,:) = tab(i,:) - tab(i,var_e)*tab(ind_s,:);
        end
    end
    
    A = tab(1:m,1:n+m+1);
    b = tab(1:m,end);
    c = -tab(end,1:n+m+1)';
    jB_aux= jB;
    zo_aux = tab(end,end);
    tab(end,end) = 0;
    %%%ENTRAMOS A FASE II%%%
    
    %%PARA NO PERDER EL TABLEAU SE COPIó AQUI LA FASE 2 PARA CORRER EL PROB AUX%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n_aux = n;
    [m,n] = size(A);
    ban = 0;
    stop = false;
    
    while (~stop)
  
        if any(tab(end,1:n) > 0)
            %disp(tab)
            [~,e] = max(tab(end,1:n)); %valor de xe, indice de xe
            ind_e = find(jN==e);
            
            %     Usando R de max descenso
            if all(tab(1:m,e) <= 0)
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
                %%T = ["Iteración número: ", num2str(iter +1)];
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
    
    %%disp(tab)
    xo_aux = tab(1:m, n+1);
    m_aux = max(m,n-m);
    xo = zeros(m_aux,1);
    B = jB;
    for i = 1:m
        if B(i)<m_aux
            xo(B(i))=xo_aux(i);
        end
    end
    

    zo = tab(end, end)+zo_aux;
    xo = xo(1:n_aux,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    
    if ban==0
        if zo<=10e-4
            if ismember(n_aux+m+1,jB)
                %%es degenerado
                
            end
            A = tab(1:m,1:n-1);
            b = tab(1:m,end);
            c_orig = c_orig(1:n-1);
            c_aux = zeros(n-1,1);
            jN = jN(jN~=n_aux+m+1);
            c_aux(jN) = -c_orig(jN) + (c_orig(jB)'*A(:,jN))';
            zo_aux = c_orig(jB)'*b;
            [xo, zo, ban, iter, B] = mSimplexFaseII_eq(A,b,c_aux,jB);
            zo = zo + zo_aux;
            xo = xo(1:n_aux);
            
   
        else
            xo = [];
            ban=-1;
            fprintf('CF vacio');
        end
    else
        fprintf('problema auxiliar no acotado');
    end

end
end
