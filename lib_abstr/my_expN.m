function v = my_expN(vmin, vmax, N)
% v = all vertices of rectangle defined by vmin and vmax
% assume, vmin/vmax not empty 
v = [];
n = length(vmin);

if (n == 1)
    v = vmin:(vmax-vmin)/(N-1):vmax;
else
    v = [];
    vr = my_expN(vmin(2:n),vmax(2:n),N);
    for j = 1:1:N
        vj = (vmin(1) + (j-1)*(vmax(1)-vmin(1))/(N-1))*ones(1,N^(n-1));
        vj = [vj; vr];
        v = [v, vj];
    end
end