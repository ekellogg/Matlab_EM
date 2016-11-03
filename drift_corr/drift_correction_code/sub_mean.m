function[normM] = sub_mean(M)
   normM = M - mean(M(:));
end