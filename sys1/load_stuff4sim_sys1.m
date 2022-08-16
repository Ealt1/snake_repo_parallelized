gpart = GridPartition({1:6,1:4});

node1.p = gpart.idx2mid(gpart.cord2idx([1;3])); 
node2.p = gpart.idx2mid(gpart.cord2idx([3;1])); 
node3.p = gpart.idx2mid(gpart.cord2idx([5;3])); 
node4.p = gpart.idx2mid(gpart.cord2idx([2;3])); 
node5.p = gpart.idx2mid(gpart.cord2idx([3;2])); 
node6.p = gpart.idx2mid(gpart.cord2idx([3;3])); 
node7.p = gpart.idx2mid(gpart.cord2idx([4;3])); 

nodes = [node1, node2, node3, node4, node5, node6, node7];

load('everything.mat');
sys = load_sys1('load');

window = [1,6,1,4];
