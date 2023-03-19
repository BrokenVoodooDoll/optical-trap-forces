function S = iscomplex(A)
    S = A .* double(imag(A) == 0);
end