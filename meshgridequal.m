function y_vec = meshgridequal(a, n)
%% made by EASHA TU RAZIA as part of Final Masters project 2021 
%% MEng Aeronautical Engineering, 2021
% a is line vector that needs to be extended
% n is the dimension

l = length(a);
y = ones( l*ones(1, n) ).*a;

y_vec = ones(l^n ,n);
b = 1:n;
for i=1:n
    b = [b(2:n), b(1)];
    temp = permute(y, b);
    y_vec(:,i) = temp(:);
end

end