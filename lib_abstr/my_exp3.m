function v = my_exp3(vmin, vmid, vmax)
% v = all vertices of rectangle defined by vmin and vmax
% assume, vmin/vmax not empty 
v = [];
n = length(vmin);

if (n == 1)
    v = [vmin, vmid, vmax];
else
    vm = vmin(1)*ones(1,3^(n-1));
    vd = vmid(1)*ones(1,3^(n-1));
    vM = vmax(1)*ones(1,3^(n-1));
    vr = my_exp3(vmin(2:n),vmid(2:n),vmax(2:n));
    vm = [vm; vr];
    vd = [vd; vr];
    vM = [vM; vr];
    v = [vm, vd, vM];
end
