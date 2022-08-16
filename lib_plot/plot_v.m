function h = plot_v(v)

if strcmp(v.type,'large')
    v_large = Polyhedron('lb', v.xy-[0.06;0.06], 'ub', v.xy+[0.06;0.06]); 
    h = plot(v_large,'color','blue', 'linestyle','-'); 
else
    v_small = Polyhedron('lb', v.xy-[0.04;0.04], 'ub', v.xy+[0.12;0.04]); 
    h1 = plot(v_small,'color',[1,1,1], 'linestyle','-'); 
    hold on;
    v_small = Polyhedron('lb', v.xy-[0.04;0.04], 'ub', v.xy+[-0.04+0.16*(max(v.f,0)/v.fmax);0.04]); 
    h2 = plot(v_small,'color',[0.5,0.5,1], 'linestyle','none'); 
    if v.f <= 0
        h3 = text(v.xy(1)+0.15, v.xy(2), ['life:', num2str(2-floor(abs(v.f)))],'Fontsize',10, 'Color', [0.5,0.5,1]);
    else
        h3 = [];
    end
    if ~v.sf
        h4 = text(v.xy(1)+0.15, v.xy(2), ['back to\newlinesafe in:', num2str(v.d)],'Fontsize',10, 'Color', [1,0,0]);
    else
        h4 = [];
    end
    if v.no==1
        h5= text(v.xy(1)-0.25, v.xy(2), 'v1' ,'Fontsize',10, 'Color', [0,0,0]);
    else
        h5= text(v.xy(1)-0.25, v.xy(2), 'v2' ,'Fontsize',10, 'Color', [0,0,0]);
    end
    h = [h1, h2, h3, h4, h5];
end




