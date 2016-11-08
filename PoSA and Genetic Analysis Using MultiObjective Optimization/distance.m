function d = distance(v1,v2)

d = [];
N1 = length(v1);
N2 = length(v2);

if N1 ~= N2
    disp('Vectors in input have different size')
    return
end

diff=zeros(N1,1);
for i=1:N1
    diff(i)=(v1(i)-v2(i))^2;
end

d=sqrt(sum(diff));