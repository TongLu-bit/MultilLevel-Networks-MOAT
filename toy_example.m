clc;clear;close all

%% Step 1 & 2
%%%% To see the initial two steps that generate the particial correlation matrix bettwen SIs and FCs, 
%%%% see "simulate_inputdata.m"
%%%% Or, one could download "W_input.mat" and use it directly here.

%% Step3: Apply MOAT to the FC-SI association matrix (500*4950)

%%%%%%%%%%%%%%%%%%%%%%%%% Step 3.1 -Detect B1(G1): %%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Step 3.1.1 Screening, to exclude most falst positive edges
ratio_threshold=99
Wp_orig=w_perm_abs;
threshold = prctile(Wp_orig(:), ratio_threshold);
Wp = Wp_orig;
Wp(Wp_orig<threshold) = 0;
% figure;imagesc(Wp);colormap jet;
% xlabel("FC",'FontSize',20,'FontWeight','bold','Color','b');
% ylabel("SI",'FontSize',20,'FontWeight','bold','Color','r');

% reorder Wp to a thresholded matrix, a smaller one.
rowsum = sum(Wp,2)'; %figure;hist(rowsum,50); title("rowsum");
colsum = sum(Wp, 1); %figure;hist(colsum,50); title("colsum");
% determine nonzero row and columns
rowsum_thrd=5.9
colsum_thrd=0.75
s_in0 = find(rowsum>rowsum_thrd);   s_out0 = 1:size(Wp_orig,1);  s_out0(s_in0) = [];  %%%SNPs
t_in0 = find(colsum>colsum_thrd);   t_out0 = 1:size(Wp_orig,2);  t_out0(t_in0) = [];  %%%fc
 
Wp_0 = Wp(s_in0, t_in0);  %342*2509
figure;imagesc(Wp_0);colorbar;
xlabel("simulated FC",'FontSize',20,'FontWeight','bold','Color','r');
ylabel("simulated SI",'FontSize',20,'FontWeight','bold','Color','r');

%%%% Step 3.1.2 Apply MOAT
r= prctile(w_perm_abs(:), ratio_threshold)
c_vec = 0.1%:0.3:5;
lambda_w=1.4;
lambda_u=1.4;
num_region=100;


tic
[s_in t_in c_max denWU_max denW_max]=greedy_multilayer(w_perm_abs,ratio_threshold, rowsum_thrd, colsum_thrd, c_vec,lambda_w,lambda_u, edge_node_list,num_region);
toc

true_t = sort(t_in0(t_in));
true_s = sort(s_in0(s_in));  %compare with true_s_idx1

N=size(Wp_0,1)
M=size(Wp_0,2)
s_out = setdiff(1:N,s_in);
t_out = setdiff(1:M,t_in);
% %figure;imagesc(w(s_in,t_in));colormap jet;colorbar;
% figure;imagesc(Wp_0([s_in sort(s_out)],[t_in sort(t_out)]));colorbar;
% xlabel("simulated FC",'FontSize',20,'FontWeight','bold','Color','r');
% ylabel("simulated SI",'FontSize',20,'FontWeight','bold','Color','r');

W_in=Wp_0([s_in sort(s_out)],[t_in sort(t_out)]);
figure;hist(reshape(W_in,1,353*2656),200);


Wp1=w_perm_abs;Wp1(w_perm_abs<0) = 0;
W_reorder=Wp1([s_in0([s_in sort(s_out)]), sort(s_out0)],[t_in0([t_in sort(t_out)]), sort(t_out0)]);

figure;imagesc(W_reorder);colorbar;colormap jet;
xlabel("simulated FC",'FontSize',20,'FontWeight','bold','Color','r');
ylabel("simulated SI",'FontSize',20,'FontWeight','bold','Color','r');


%%%% set extracted node-pair =1
extracted_i=edge_node_list(true_t,2);
extracted_j=edge_node_list(true_t,3);
see=[extracted_i extracted_j];


u1=zeros(100,100);
for i=1:size(extracted_i,1)
    u1(extracted_i(i),extracted_j(i))=2;
end 

u1=u1+u1.';
figure;imagesc(u1);ax=gca;ax.FontSize=18;
xlabel('Brain Regions','FontSize',24);
ylabel('Brain Regions','FontSize',24);

%density of B1(G1)
denW_max 
%density of B(G)
den_raw_w=length(find(w_perm_abs>r))/(sqrt(size(w_perm_abs,1)*size(w_perm_abs,2)))^lambda_w


%%% Since density of B1(G1) > density of B1(G1), we now proceed to detect
%%% 2nd dense sub-network, B2(G2)



%%%%%%%%%%%%%%%%%%%%%%%%% Step 3.2 -Detect B2(G2): %%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Step 3.2.1 fill the medium of B(G) into B1(G1)
med=median(w_perm_abs,'all')
w_fill=w_perm_abs;
w_fill(true_s, true_t)=med;
figure;imagesc(w_fill);ax=gca;ax.FontSize=18;
xlabel("FC");ylabel("SI");


%%%% Step 3.2.2 Screening, to exclude most falst positive edges
threshold = prctile(w_fill(:), ratio_threshold)
Wp = w_fill;
Wp(w_fill<threshold) = 0;

rowsum = sum(Wp,2)'; figure;hist(rowsum,50); title("rowsum");
colsum = sum(Wp, 1); figure;hist(colsum,50); title("colsum");
% determine nonzero row and columns
rowsum_thrd2=10
colsum_thrd2=1.65
s_in02 = find(rowsum>rowsum_thrd2);   s_out02 = 1:size(w_fill,1);  s_out02(s_in02) = [];  %%%SNPs
t_in02 = find(colsum>colsum_thrd2);   t_out02 = 1:size(w_fill,2);  t_out02(t_in02) = [];  %%%fc
Wp_02 = Wp(s_in02, t_in02);  %342*2509 

c_vec = 0.1%:1:5;
lambda_w=1.4;
lambda_u=1.4;
num_region=100;


tic
[s_in2 t_in2 c_max2 denWU_max2 denW_max2]=greedy_multilayer(w_fill,ratio_threshold, rowsum_thrd2, colsum_thrd2, c_vec,lambda_w,lambda_u, edge_node_list,num_region);
toc


true_t2 = sort(t_in02(t_in2));
true_s2 = sort(s_in02(s_in2));  %compare with true_s_idx2

N=size(Wp_02,1)
M=size(Wp_02,2)
s_out2 = setdiff(1:N,s_in2);
t_out2 = setdiff(1:M,t_in2);
%figure;imagesc(w_fill(s_in1,t_in1));colormap jet;colorbar;
% figure;imagesc(Wp_02([s_in2 sort(s_out2)],[t_in2 sort(t_out2)]));colorbar;
% xlabel("FC",'FontSize',16,'FontWeight','bold','Color','r');
% ylabel("SI",'FontSize',16,'FontWeight','bold','Color','r');

Wp2=w_fill;Wp2(w_fill<0) = 0;
W_reorder2=Wp2([s_in02([s_in2 sort(s_out2)]), sort(s_out02)],[t_in02([t_in2 sort(t_out2)]), sort(t_out02)]);

figure;imagesc(W_reorder2);colorbar;colormap jet;
xlabel("simulated FC",'FontSize',20,'FontWeight','bold','Color','r');
ylabel("simulated SI",'FontSize',20,'FontWeight','bold','Color','r');


%%%% set extracted node-pair =1
extracted_i2=edge_node_list(true_t2,2);
extracted_j2=edge_node_list(true_t2,3);
see=[extracted_i2 extracted_j2];


u2=zeros(100,100);
for i=1:size(extracted_i2,1)
    u2(extracted_i2(i),extracted_j2(i))=1;
end 

u2=u2+u2.';
figure;imagesc(u2);ax=gca;ax.FontSize=18;
xlabel('Brain Regions','FontSize',24);
ylabel('Brain Regions','FontSize',24);
%[cindx,comp_size]=conncomp(graph(w))



%% Step4: Detected Results Visulization 
s_rest=setdiff(1:500,[true_s true_s2]);
t_rest=setdiff(1:4950,[true_t true_t2]);

figure;imagesc(w_perm_abs([true_s true_s2 sort(s_rest)],[true_t true_t2 sort(t_rest)]));colorbar;colormap jet;
ax=gca;ax.FontSize=18;ax.FontWeight='bold';
xlabel("FC",'FontSize',20,'FontWeight','bold','Color','k');
ylabel("SI",'FontSize',20,'FontWeight','bold','Color','k');

u_all=u1+u2;
figure;imagesc(u_all);ax=gca;ax.FontSize=18;ax.FontWeight='bold';colorbar;
xlabel('Brain Regions','FontSize',20,'FontWeight','bold','Color','k');
ylabel('Brain Regions','FontSize',20,'FontWeight','bold','Color','k');
