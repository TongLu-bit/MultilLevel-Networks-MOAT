%% Step 1: Simulate SI, FC measures respectively 
%%%% 1. 500 SI measures and 4950 FC measures will be simulated, where FCs
%%%% compose a brain network with 100 regions
%%%% 2. Two dense subnetworks B1(G1), B2(G2) with remarkbaly high FCN-SI
%%%% partial corr. are generated within B(G)

%%%% Define effect sizes 
b1=0.7;b2=0.4;
noise1=0.2;noise2=3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 1.1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Construct a (500+4950)*(500+4950) corr. matrix for an MVN where SI, FC
%%%% can be simulated from
N=500+4950;
N1=N*(N-1)/2;
r = noise1.*rand(N1,1);  %scalar, noise1*uniform(0,1)
sig=squareform(r);
figure;imagesc(sig);ax=gca;ax.FontSize=18;


%%%% Subnetwork 1: B1(G1)
%%%% B1 consists of 40 SI measures and 435 FC measures that compose G1 of 30 brain regions
%%%% Therefore, there are 40+435=475 variables to be considered. 
length1=475
t1=zeros(length1,length1);
for i=1:length1
    for j=1:length1
        if j>i
            t1(i,j)=0.3+0.6*unifrnd(b1,1);  %uniform(b1,1)
        end 
    end 
end 
%sanity check : if all off-diagonal elements <=1
find(t1>=1);

d=ones(1,length1)';
block1 = diag(d)+t1+t1.';
%figure; imagesc(block1);ax=gca;ax.FontSize=18;
sig(1:length1,1:length1)=block1;
figure;imagesc(sig);ax=gca;ax.FontSize=18;

%%%% Subnetwork 2: B2(G2)
%%%% B2 consists of 60 SI measures and 190 FC measures that compose G2 of 20 brain regions
%%%% Therefore, there are 60+190=250 variables to be considered. 
length2=250;
t2=zeros(length2,length2);
for i=1:length2
    for j=1:length2
        if j>i
            t2(i,j)=0.3+0.6*unifrnd(b2,1);  %uniform(b2,1)
        end 
    end 
end 
%sanity check : if all off-diagonal elements <=1
find(t2>=1);

d=ones(1,length2)';
block2 = diag(d)+t2+t2.';
%figure; imagesc(block2);ax=gca;ax.FontSize=18;
sig((length1+1):(length1+length2),(length1+1):(length1+length2))=block2;
figure;imagesc(sig);ax=gca;ax.FontSize=18;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 1.2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Add extra noise, if prefer. 
% A = randi(noise2,1275,1275);  %returns a 1275*1275 matrx of integres between 1 and 10
% A=A/10;  %constrain the range to (0.1,0.5)
% max(max(A))
% 
% B = randi(noise2,1275,1275);  %returns a 1275*1275 matrx of integres between 1 and 10
% B=-1*B/10;  %constrain the range to (0.1,0.5)
% max(max(B))
% 
% C=A+B;
% %C(randsample(320,250),randsample(320,250))=0;
% 
% C(1:320,1:320)=0;
% max(max(C))
% %figure;imagesc(C);ax=gca;ax.FontSize=18;
% 
% %insert 0 to random (i,j) position, control the noise ratio:
% num_non_noise_voxel=1000  %about 95% of 1275
% noise_i=randsample(1275, num_non_noise_voxel);
% noise_j=randsample(1275, num_non_noise_voxel);
% C(noise_i,noise_j)=0;
% %figure;imagesc(C);ax=gca;ax.FontSize=18;
% 
% sig_noist=sig+C;
% figure;imagesc(sig_noist);ax=gca;ax.FontSize=18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 1.3: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% solve SPD issue of the sig matrix, if any 
%sig1=nearestSPD(sig_noist);

sig1=nearestSPD(sig);
figure;imagesc(sig1);ax=gca;ax.FontSize=18;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 1.4: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Now, generate 500 SIs, 4950 FCs from constructed MVN

mu=zeros(1,N)'; % Mean of MVN 
dat = mvnrnd(mu, sig1, 100); %100*5450, generate 100 subjects

%%%% Obtain the corresponding 500 SIs from dat 
si=[dat(:,1:40),dat(:,476:535),dat(:,726:1125)]; %sub * value

%%%% Obtain the corresponding 4950 FCs from dat 
fc1=dat(:,41:475);  %FCs in B1(G1)
fc2=dat(:,536:725); %FCs in B2(G2)
fc_rest=dat(:,1126:end);
fc_all=[fc1 fc2 fc_rest];


%% Step2: Simulate partial corr. matrix for FC-SI associations (500*4950)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 2.1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Ground truth
w_true=corr(si,fc_all);
w_true0=abs(w_true);
figure;imagesc(w_true0); colormap jet;colorbar;ax=gca;ax.FontSize=18;ax.FontWeight='bold';
xlabel("FC",'FontSize',20,'FontWeight','bold','Color','k');
ylabel("SI",'FontSize',20,'FontWeight','bold','Color','k');

figure;hist(reshape(w_true0,1,500*4950),200);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 2.2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Inner layer network structure visulization, G
g=zeros(100,100);
g(1:30,1:30)=2;
g(31:50,31:50)=1;
figure;imagesc(g);colorbar;ax=gca;ax.FontSize=18;ax.FontWeight='bold';
xlabel("Brain Regions",'FontSize',20,'FontWeight','bold','Color','k');
ylabel("Brain Regions",'FontSize',20,'FontWeight','bold','Color','k');


%%%% Construct n*n inner layer network, G:
edge_node_list=zeros(4950,3);
edge_node_list(:,1)=1:4950;

x=[];y=[];
for i=1:100
    for j=1:100        
        if j>i
            x=[x i];
            y=[y j];
        end 
    end 
end 
edge_node_list(:,2)=x';
edge_node_list(:,3)=y';

%edge index of voxels in u1 (the dense sub-area in U)
u1_idx=find(edge_node_list(:,2)>=1 & edge_node_list(:,2)<=30 & edge_node_list(:,3)>=1 & edge_node_list(:,3)<=30);
u2_idx=find(edge_node_list(:,2)>=31 & edge_node_list(:,2)<=50 & edge_node_list(:,3)>=31 & edge_node_list(:,3)<=50);


%assign that first abnormal 190 columns in fc to u1_idx:
fc_perm=zeros(100,4950);
fc_perm(:,u1_idx)=fc1;
fc_perm(:,u2_idx)=fc2;


neg_idx=find(sum(fc_perm)==0);
fc_perm(:,neg_idx)=fc_rest;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 2.2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 500*4950 partial corr. matrix (FC nodes permuted, SI nodes not yet)
w=corr(si,fc_perm);%correlation between **columns** of X and Y.
figure;imagesc(w);ax=gca;ax.FontSize=18;
xlabel("FC");ylabel("SI");

w_abs=abs(w);

%%%% Permute SIs
row_idx=randperm(size(w,1));
true_s_idx1=[];
for i=1:40
    true_s_idx1=[true_s_idx1 find(row_idx==i)];  %for reference, to compare with true_s1 yileded by greedy later
end
true_s_idx1=sort(true_s_idx1);

true_s_idx2=[];
for i=41:100
    true_s_idx2=[true_s_idx2 find(row_idx==i)];  %for reference, to compare with true_s1 yileded by greedy later
end
true_s_idx2=sort(true_s_idx2);


w_perm=w(row_idx,:);
figure;imagesc(w_perm);colorbar;ax=gca;ax.FontSize=18;
xlabel("FC");ylabel("SI");title("input")

%%%% Finally, input data for MOAT
w_perm_abs=abs(w_perm);
figure;imagesc(w_perm_abs);colorbar;colormap jet;colorbar;ax=gca;ax.FontSize=18;ax.FontWeight='bold';
xlabel("FC",'FontSize',20,'FontWeight','bold','Color','k');
ylabel("SI",'FontSize',20,'FontWeight','bold','Color','k');


