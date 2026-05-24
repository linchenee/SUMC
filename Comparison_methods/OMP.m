function [x] = OMP(y,A,hatA,r)
% ----------------------------------------   
% Orthogonal matching pursuit (OMP) method
% ----------------------------------------
[~,coef] = cal_coef(y,hatA,r);
x = A*coef;

end

function [Ids,coef] = cal_coef(y,A,t)
% ---------------------------------------------------   
% This function is used to calculate the coefficient.
% ---------------------------------------------------
[y_rows,y_columns] = size(y);
if y_rows < y_columns
    y = y';
end
[M,N] = size(A);
coef = zeros(N,1);
At = zeros(M,t);
Ids = zeros(1,t);
r_n = y;
for i = 1:t
    product = A'*r_n;
    [~,id] = max(abs(product));
    At(:,i) = A(:,id);
    Ids(i) = id;
    A(:,id) = zeros(M,1);
    coef_now = (At(:,1:i)'*At(:,1:i))\(At(:,1:i)'*y);
    w = warning('query','last');
    warning('off',w.identifier);
    r_n = y - At(:,1:i)*coef_now;        
end
coef(Ids) = coef_now;

end