function vector_y=chebyshev_basis_y(n,space_y)
vector_y=[]
syms y
if n==1
    vector_y=[1]
elseif n==2
    vector_y=[1,y]
else n>2
    vector_y=[1,y]
    for i=3:n
    vector_y=[vector_y 2*y*vector_y(1,i-1)-vector_y(1,i-2)]
    end
end
end