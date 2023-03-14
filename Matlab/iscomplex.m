function S = iscomplex(A)
[n,m] = size(A);
S = zeros(n,m);

for ii = 1:n
    for jj = 1:m
        S(ii,jj) = A(ii,jj)*double(isreal(A(ii,jj)));
    end
end

end