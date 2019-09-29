function vector_x=chebyshev_basis_x(n,space_x)
vector_x=[]
syms x
if n==1
    vector_x=[1]
elseif n==2
    vector_x=[1,x]
else n>2
    vector_x=[1,x]
    for i=3:n
    vector_x=[vector_x 2*x*vector_x(1,i-1)-vector_x(1,i-2)]
    end
end
end