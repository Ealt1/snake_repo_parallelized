function [h0, h1, h2, ho, hol] = plot_x(x, nodes)

p0 = x(1);
v0.type = 'large';
v0.xy = nodes(p0).p + [0.1;0.35];
h0 = plot_v(v0); 

p1 = x(2);
p2 = x(3);
f1 = x(4);
f2 = x(5);
d1 = x(8);
d2 = x(9);

v1.type = 'small';
v1.no = 1;
v1.xy = nodes(p1).p + [0.1;0.1];
v1.f = f1;
v1.fmax = 3;
v1.sf = 1;
if ismember(p1,[1 2 3])
    v1.sf = 0;
end
v1.d = d1; 

v2.type = 'small';
v2.no = 2;
v2.xy = nodes(p2).p + [0.1;0.2];
v2.f = f2;
v2.fmax = 3;
v2.sf = 1;
if ismember(p2,[1 2 3])
    v2.sf = 0;
end
v2.d = d2;


h1 = plot_v(v1);
h2 = plot_v(v2); 

po = x(6); 
to = x(7);
if ~(po == p1 || po == p2)
    ho = scatter(nodes(po).p(1), nodes(po).p(2), 3, [1, 0, 0], 'Linewidth', 5);
    hol = text(nodes(po).p(1), nodes(po).p(2)+0.25, ['time left:', num2str(to)],'Fontsize',12, 'Color', [1 0 0]);
else
    ho = scatter(nodes(po).p(1), nodes(po).p(2), 3, [0.8, 0.8, 0.8], 'Linewidth', 5);
    hol = [];
end



