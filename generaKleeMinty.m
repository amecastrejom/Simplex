function[c, A, b] = generaKleeMinty(m)
 c = (-1)*ones(m,1);
 A = eye(m) + 2*tril(ones(m),-1);
 b = ones(m,1);
 for j = 1:m
     b(j) = 2^j - 1;
 end
 b;

end