function v = my_exp(vmin, vmax)
% v = all vertices of rectangle defined by vmin and vmax
% assume, vmin/vmax not empty 
v = [];
n = length(vmin);

if (n == 1)
    v = [vmin, vmax];
else
    vm = vmin(1)*ones(1,2^(n-1));
    vM = vmax(1)*ones(1,2^(n-1));
    vr = my_exp(vmin(2:n), vmax(2:n));
    vm = [vm; vr];
    vM = [vM; vr];
    v = [vm, vM];
end
