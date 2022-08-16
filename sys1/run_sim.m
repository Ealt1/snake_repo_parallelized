sys = load_sys1();
load('everything.mat');

addpath('plots')
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
end

axis equal
axis([1,6,1,4])

% Run sim
xt_idx = C(1); 
while 1
    if ~ismember(xt_idx, C)
        error('not in cinv!?');
        break;
    end
    xt = sys.gpart_x.idx2mid(xt_idx);
    [h0, h1, h2, ho] = plot_x(xt, nodes);
    axis([1,6,1,4])
    % prompt = 'continue?';
    % not_important = input(prompt); 
    % make decision
    % u = K{find(C == xt_idx)}(1);
    u = least_move(xt, K{find(C == xt_idx)}, sys.gpart_u);
    xt_idxs_plus = find(sys.tau{u}(xt_idx,:)); 
    if length(xt_idxs_plus) >= 2
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
    delete([h0, h1, h2]);
    xt_plus = sys.gpart_x.idx2mid(xt_idx);
    plot_trans(xt, xt_plus, nodes); 
    delete([ho]);
end

