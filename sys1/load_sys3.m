function sys = load_sys3(flag)

if nargin < 1
    flag = '';
end

n = 7;
Map = zeros(n,n);
Map(1,4) = 1;
Map(4,5) = 1;
Map(4,6) = 1;
Map(5,2) = 1;
Map(5,6) = 1;
Map(6,7) = 1;
Map(7,3) = 1;
Map = (Map + Map' + eye(n)) > 0.5; 

sys.M = Map;

sys.p_large_safe = [4 5 6 7];
sys.p_target = [1 2 3];
l = length(sys.p_target);

sys.f_max = 3;
sys.d_max = 2;
sys.to_life = 2; 

% x = [p_large, p_1, p_2, f_1, f_2, p_o, t_o, d_1, d_2]
grid_x = {1-0.5:n+0.5,...
          1-0.5:n+0.5,...
          1-0.5:n+0.5,...
          -2-0.5:sys.f_max+0.5,...
          -2-0.5:sys.f_max+0.5,...
          1-0.5:l+0.5,...
          0-0.5:sys.to_life+0.5,...
          0-0.5:sys.d_max+0.5,...
          0-0.5:sys.d_max+0.5}; 
sys.gpart_x = GridPartition(grid_x);

grid_u = {1-0.5:n+0.5,...
          1-0.5:n+0.5,...
          1-0.5:n+0.5}; 
sys.gpart_u = GridPartition(grid_u);


% load transitions
if strcmp(flag, 'load')
load('sys_tau_1to100.mat', 'sys_tau_1to100');
load('sys_tau_101to200.mat', 'sys_tau_101to200');
load('sys_tau_201toend.mat', 'sys_tau_201toend');
sys.tau = [sys_tau_1to100, sys_tau_101to200, sys_tau_201toend];
return; 
end


N = length(sys.gpart_x) + 1
M = length(sys.gpart_u);

% compute transitions
% tau = {};
% for k = 1:1:M
%     tau_k = sparse(N,N);
%     tau_k(N, N) = 1;
%     tau{end+1} = tau_k;
% end

tau = sptensor([],[],[N N M]);
tau(N,N,:) = 1;

progressStepSize = 1000;
ppm = ParforProgMon('Example', N, progressStepSize);
my_i_k = zeros(N,1000);


test = zeros(N,1);
parfor i = 1:1:N-1
%parfor i = 1:1:H-1
%     if ~mod(i,1000)
%         i
%     end
    if mod(i,progressStepSize)==0
        ppm.increment();
    end
    x_i = sys.gpart_x.idx2mid(i);
    p0 = x_i(1);
    p1 = x_i(2);
    p2 = x_i(3);
    f1 = x_i(4);
    f2 = x_i(5);
    po = x_i(6); 
    to = x_i(7);
    d1 = x_i(8);
    d2 = x_i(9);
    
    P0_plus = find(Map(p0,:));
    P1_plus = find(Map(p1,:));
    P2_plus = find(Map(p2,:));
   
    mytemp = zeros(1,1000);
    for p0_plus = P0_plus
        for p1_plus = P1_plus
            for p2_plus = P2_plus
                k = sys.gpart_u.pnt2idx([p0_plus;p1_plus;p2_plus]);
                
                unsafe_ik = 0;
                if ~ismember(p0, sys.p_large_safe)    
                    mytemp(k) = 1;    
                    %tau{k} = my_val;
                    
                else
                    if Map(p0,p0_plus) && Map(p1,p1_plus) && Map(p2,p2_plus) 
                        % valid move on the Map
                        
                        % dead time update
                        if ismember(p1,setdiff(1:n,sys.p_large_safe))
                            d1_plus = d1 - 1;
                        else
                            d1_plus = sys.d_max;
                        end
                        if ismember(p2,setdiff(1:n,sys.p_large_safe))
                            d2_plus = d2 - 1;
                        else
                            d2_plus = sys.d_max;
                        end
                        if d1_plus <= -0.5 || d2_plus <= -0.5
                            mytemp(k) = 1;
                            %tau{k}(i,N) = 1; 
                            unsafe_ik = 1;
                        end

                        % fuel update
                        if ~unsafe_ik
                            if f1 >= 0.5 % enough fuel for nontrivial move
                                if p1 ~= p1_plus
                                    f1_plus = f1 - 1;
                                else
                                    f1_plus = f1;
                                end
                            else % zero fuel
                                if p1 ~= p1_plus
                                    mytemp(k) = 1;
                                    %tau{k}(i,N) = 1; 
                                    unsafe_ik = 1;
                                else
                                    f1_plus = f1 - 1;
                                end
                            end

                            if f2 >= 0.5
                                if p2 ~= p2_plus  
                                    f2_plus = f2 - 1;
                                else
                                    f2_plus = f2;
                                end
                            else
                                if p2 ~= p2_plus  
                                    mytemp(k) = 1;
                                    %tau{k}(i,N) = 1;
                                    unsafe_ik = 1;
                                else
                                    f2_plus = f2 - 1;
                                end
                            end
                        end

                        if ~unsafe_ik
                            if f1_plus <= -2.5 || f2_plus <= -2.5
                                % over-draw fuel
                                mytemp(k) = 1;
                                %tau{k}(i, N) = 1;
                                % unsafe_ik = 1;
                            else
                                % re-fuel from the large vehicle "v0"
                                if p1_plus == p0_plus
                                    f1_plus = sys.f_max;
                                end
                                if p2_plus == p0_plus
                                    f2_plus = sys.f_max;
                                end

                                if to >= -0.5 && (p1 == po || p2 == po)
                                    % capture now
                                    to_plus = sys.to_life; 
                                    po_plus = sys.p_target;
                                else
                                    % capture later
                                    to_plus = to - 1;
                                    po_plus = po;
                                end

                                if to_plus <= -0.5
                                    % target dies before being captured
                                    mytemp(k) = 1;
                                    %tau{k}(i, N) = 1;
                                    % unsafe_ik = 1;
                                else
                                    for p = po_plus
                                        i_plus = sys.gpart_x.pnt2idx([p0_plus; p1_plus; p2_plus; f1_plus; f2_plus; p; to_plus; d1_plus; d2_plus])    
                                    end
                                    fprintf('size(p) is %s\n', mat2str(size(po_plus)))
                                    %test(i) = 1;
                                    %tau{k}(i, i_plus) = 1;
                                end
                            end
                        end
                    else
                        % invalid move on the Map
                        mytemp(k) = 1;
                        %tau{k}(i, N) = 1;
                    end
                end
            end
        end
    end
my_i_k(i,:) = mytemp;    
end
ppm.delete();


%ind = find(my_i_k==1);
% parfor my_ind = ind
%     [i, k] = ind2sub(size(my_i_k),my_ind);
%     tau{k}(i,N) = 1;
% end
% 
% find(i
%tau_mat = cell2mat(permute(tau,[1,3,2]));
% tau_mat =  cell2mat(reshape(tau,1,1,[]));
% 
%mytest = size(find(test==1))
ind = find(my_i_k==1);
[i, k] = ind2sub(size(my_i_k),ind);
myN = N*ones(size(i));
size(N)
size(i)
tau([i,myN,k]) = 1;

% for i = 1:N
% %     if ~mod(i,1000)
% %         i
% %     end
%     k = find(my_i_k(i,:)==1);
% %     tau_temp = tau(i,:,:);
% %     tau_temp(i,N) = 1;
%         tau(i,N,k) = 1;
% end

% ind = find(my_i_k==1);
% [i, k] = ind2sub(size(my_i_k),ind);

% for h = 1:length(i)
%     tau{k(h)}(i(h),N) = 1;
% end



% for i = 1:N
%     k = find(my_i_k(i,:)==1);
%     for myk = k
%         tau{myk}(i,N) = 1;
%     end
% end
sys.tau = tau;
