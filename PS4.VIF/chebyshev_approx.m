function cheby=chebyshev_approx(x,n)
X=x(:);
lx=size(X,1);

    if n<0;
        error('n should be a positive integer');
    end
    
switch n;

    case 0;
        cheby=[ones(lx,1)];
    case 1;
        cheby=[ones(lx,1) X];
    otherwise
        cheby=[ones(lx,1) X];
    for i=3:n+1;
        cheby=[cheby 2*X.*cheby(:,i-1)-cheby(:,i-2)];
    end
end