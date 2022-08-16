function plot_trans(x,x_plus, nodes, window, N)

if nargin < 4
    window = [1,6,1,4]; 
end

if nargin < 5
    N = 12;
end

p0 = x(1);
p0_plus = x_plus(1);
pp0 = linspace(nodes(p0).p(1), nodes(p0_plus).p(1),N+2) + 0.1;
pp0 = [pp0; linspace(nodes(p0).p(2), nodes(p0_plus).p(2),N+2) + 0.35];
pp0 = pp0(:,2:end-1);

p1 = x(2);
p1_plus = x_plus(2);
pp1 = linspace(nodes(p1).p(1), nodes(p1_plus).p(1),N+2) + 0.1;
pp1 = [pp1; linspace(nodes(p1).p(2), nodes(p1_plus).p(2),N+2) + 0.1];
pp1 = pp1(:,2:end-1);


p2 = x(3);
p2_plus = x_plus(3);
pp2 = linspace(nodes(p2).p(1), nodes(p2_plus).p(1),N+2) + 0.1;
pp2 = [pp2; linspace(nodes(p2).p(2), nodes(p2_plus).p(2),N+2) + 0.2];
pp2 = pp2(:,2:end-1);

f1 = x(4);
if p0_plus == p1_plus
    f1_plus = f1; 
else
    f1_plus = x_plus(4);
end
ff1 = linspace(f1, f1_plus,N+2);
ff1 = ff1(2:end-1);

f2 = x(5);
if p0_plus == p2_plus
    f2_plus = f2; 
else
    f2_plus = x_plus(5);
end
ff2 = linspace(f2, f2_plus,N+2);
ff2 = ff2(2:end-1);

po = x_plus(6); 

d1 = x(8);
d2 = x(9);

for i = 1:1:N
    v0.type = 'large';
    v0.xy = pp0(:,i);
    h0 = plot_v(v0);
    
    v1.type = 'small';
    v1.no = 1;
    v1.xy = pp1(:,i);
    v1.f = ff1(i);
    v1.fmax = 3;
    v1.sf = 1;
    if ismember(p1,[1 2 3])
        v1.sf = 0;
    end
    v1.d = d1;
    h1 = plot_v(v1);
    
    v2.type = 'small';
    v2.no = 2;
    v2.xy = pp2(:,i);
    v2.f = ff2(i);
    v2.fmax = 3;
    v2.sf = 1;
    if ismember(p2,[1 2 3])
        v2.sf = 0;
    end
    v2.d = d2;
    h2 = plot_v(v2);
    
    if i<=N/2
        ho = scatter(nodes(po).p(1), nodes(po).p(2), 0.1, [0.8, 0.8, 0.8], 'Linewidth', 0.1);
    else
        ho = scatter(nodes(po).p(1), nodes(po).p(2), 3, [1, 0, 0], 'Linewidth', 5*(i-N/2)/(N/2));
    end
    axis(window)
    % pause(0.01)
    % pause(0.1/N)
    pause(1e-5)
    delete([h0, h1, h2, ho]);
end


