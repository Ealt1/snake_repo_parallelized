function v = my_expNs(vmin, vmax, Ns)
% v = all vertices of rectangle defined by vmin and vmax
% assume, vmin/vmax not empty 
v = [];
n = length(vmin);

Ns = Ns + (Ns==1);
if (n == 1)
    if (vmin ~= vmax)
        v = vmin:(vmax-vmin)/(Ns(1)-1):vmax;
    else
        v = [vmin vmax];
    end
else
    v = [];
    vr = my_expNs(vmin(2:n),vmax(2:n),Ns(2:n));
    for j = 1:1:Ns(1)
        M = prod(Ns(2:n));
        vj = (vmin(1) + (j-1)*(vmax(1)-vmin(1))/(Ns(1)-1))*ones(1,M);
        vj = [vj; vr];
        v = [v, vj];
    end
end




