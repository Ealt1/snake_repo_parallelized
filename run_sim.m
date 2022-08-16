clear,clc

addpath('lib_abstr')
addpath('lib_cinv_synthesis')
addpath('lib_plot')
addpath('sys1')

load('sys1/everything.mat'); 
sys = load_sys1('load');

gpart = GridPartition({1:6,1:4});

node1.p = gpart.idx2mid(gpart.cord2idx([1;3])); 
node2.p = gpart.idx2mid(gpart.cord2idx([3;1])); 
node3.p = gpart.idx2mid(gpart.cord2idx([5;3])); 
node4.p = gpart.idx2mid(gpart.cord2idx([2;3])); 
node5.p = gpart.idx2mid(gpart.cord2idx([3;2])); 
node6.p = gpart.idx2mid(gpart.cord2idx([3;3])); 
node7.p = gpart.idx2mid(gpart.cord2idx([4;3])); 

nodes = [node1, node2, node3, node4, node5, node6, node7];

% Plot Map
hold on
for i = 1:1:length(nodes)
    scatter(nodes(i).p(1), nodes(i).p(2), 3, [0.8, 0.8, 0.8], 'Linewidth', 10);
    for j = find(sys.M(i,:))
        plot([nodes(i).p(1), nodes(j).p(1)], [nodes(i).p(2), nodes(j).p(2)], 'Linewidth', 1,'Color',[0.8, 0.8, 0.8]);
    end
    
    if ismember(i, [1,2,3])
        plotCell_by_idx(gpart, gpart.pnt2idx(nodes(i).p), 1, 2, [1 0 0])
    end
end

for i = 1:1:length(nodes)
    addTxt(gpart, gpart.pnt2idx(nodes(i).p), num2str(i), 'm', [0.7, 0.7, 0.7], 12); 
end

axis equal
grid off 
set(gca,'xtick',[])
set(gca,'ytick',[])
axis([1,6,1,4])

% Run sim
appearing = 0; 
xt_idx = C(1); 
while 1
    if ~ismember(xt_idx, C)
        error('not in cinv!?');
        break;
    end
    xt = sys.gpart_x.idx2mid(xt_idx);
    [h0, h1, h2, ho, hol] = plot_x(xt, nodes);
    axis([1,6,1,4])
    % prompt = 'continue?';
    % not_important = input(prompt); 
    % make decision
    % u = K{find(C == xt_idx)}(1);
    u = least_move(xt, K{find(C == xt_idx)}, sys.gpart_u);
    xt_idxs_plus = find(sys.tau{u}(xt_idx,:)); 
    
    h_info_c = [];
    if (xt(2) == xt(6)) || (xt(3) == xt(6))
        h_info_c = text(1, 4.5, 'item captured','Fontsize',15, 'Color', [0,0,0]);
    end
    h_info_ch = [];
    if xt(1) == xt(2) || xt(1) == xt(3)
        if xt(1) ~= xt(3)
            h_info_ch = text(1, 4.2, 'small vehicle v1 recharged','Fontsize',15, 'Color', [0,0,0]);
        else
            if xt(1) ~= xt(2)
                h_info_ch = text(1, 4.2, 'small vehicle v2 recharged','Fontsize',15, 'Color', [0,0,0]);
            else
                h_info_ch = text(1, 4.2, 'small vehicles v1 v2 recharged','Fontsize',15, 'Color', [0,0,0]);
            end
        end
    end
    
    
    if length(xt_idxs_plus) >= 2
        appearing = 1;
        prompt = 'choose next item {1 2 3}: ';
        po_plus = input(prompt); 
        for i = xt_idxs_plus
            xt_potential = sys.gpart_x.idx2mid(i);
            if xt_potential(6) == po_plus
                xt_idx = i;
                break;
            end
        end
    else
        xt_idx = xt_idxs_plus(1); 
    end
    xt_plus = sys.gpart_x.idx2mid(xt_idx);
    
    % show notifications:
    if appearing
        h_info_o = text(1, 4.5, 'new item appearing','Fontsize',15, 'Color', [0,0,0]);
        appearing = 0;
    else
        if xt_plus(6)==1
             h_info_o = text(1, 4.5, 'item appears at node 1','Fontsize',15, 'Color', [0,0,0]);
        end
        if xt_plus(6)==2
             h_info_o = text(1, 4.5, 'item appears at node 2','Fontsize',15, 'Color', [0,0,0]);
        end
        if xt_plus(6)==3
             h_info_o = text(1, 4.5, 'item appears at node 3','Fontsize',15, 'Color', [0,0,0]);
        end
    end
    
    delete([h_info_c, h_info_ch]);
    h_info_ch = [];
    if xt(1) == xt(2) || xt(1) == xt(3)
        if xt(1) ~= xt(3)
            h_info_ch = text(1, 4.2, 'small vehicle v1 recharged','Fontsize',15, 'Color', [0,0,0]);
        else
            if xt(1) ~= xt(2)
                h_info_ch = text(1, 4.2, 'small vehicle v2 recharged','Fontsize',15, 'Color', [0,0,0]);
            else
                h_info_ch = text(1, 4.2, 'small vehicles v1 v2 recharged','Fontsize',15, 'Color', [0,0,0]);
            end
        end
    end
    pause(0.5)
    delete([h0, h1, h2]);
    delete([h_info_ch]);
    
    plot_trans(xt, xt_plus, nodes); 
    delete([ho]);
    delete([hol]);
    delete([h_info_o]);
end

