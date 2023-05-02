function [s_in t_in c_max denWU_max denW_max]=greedy_multilayer(Wp_orig,ratio_threshold, rowsum_thrd, colsum_thrd,c_vec,lambda_w,lambda_u, edge_node_list,num_region)
%%% ratio_threshold: percentile thershold 

    threshold = prctile(Wp_orig(:), ratio_threshold);
    Wp = Wp_orig;
    Wp(Wp_orig<threshold) = 0;
%     figure;imagesc(Wp);colormap jet;
%     xlabel("FC",'FontSize',20,'FontWeight','bold','Color','b');
%     ylabel("SNPs",'FontSize',20,'FontWeight','bold','Color','r');

    % reorder Wp to a thresholded matrix, a smaller one.
    rowsum = sum(Wp,2)';
    colsum = sum(Wp, 1);
    % determine nonzero row and columns
    s_in0 = find(rowsum>rowsum_thrd);   s_out0 = 1:size(Wp_orig,1);  s_out0(s_in0) = [];  %%%SNPs
    t_in0 = find(colsum>colsum_thrd);   t_out0 = 1:size(Wp_orig,2);  t_out0(t_in0) = [];  %%%fc
     
    Wp_0 = Wp(s_in0, t_in0);  %734*4199
    
%     figure;imagesc(Wp_0); colormap jet;
%     xlabel("FC",'FontSize',20,'FontWeight','bold','Color','b');
%     ylabel("SNPs",'FontSize',20,'FontWeight','bold','Color','r');
%     
    W=Wp_0;
    W1 = W;
    dST_vec = zeros(size(c_vec,2),2);  %1st col: stores den(W)+den(U); 2nd col: stores den(W)
    for l=1:size(c_vec,2)
        
        text = ['Current c_vec: c=',num2str(l)];
        disp(text)


        W=abs(W1);
        c=c_vec(l);
        density_list = zeros((size(W,1)+size(W,2))-1,11);
        retain_edge_list = cell((size(W,1)+size(W,2))-1, 1);

        for i=1:(size(W,1)+size(W,2))

            if rem(i,50)==0
             text1 = ['Current iteration: i=',num2str(i)];
             disp(text1)
            end 
          
            if W==0
                break
            end
            % mean of columns
            C=sum(W,1); %1*1225
            % mean of rows
            R=sum(W,2); %1*50
            %[dT,IndT]=min(C(C>0));
            % min sum in column 
            dT=min(C(C>0));
            IndT=find(C==dT);
            % min sum in row 
            dS=min(R(R>0));
            IndS=find(R==dS);
            %c=|S|/|T|=30/50
            %c=1;%sum(R~=0)/sum(C~=0);
            if sqrt(c)*dS <= dT/sqrt(c)
                W(IndS(1),:)=zeros(1,size(W,2));
                density_list(i,1)=sum(sum(W,2)~=0);
                density_list(i,2)=sum(sum(W,1)~=0);
                density_list(i,3)=1;
            else
                W(:,IndT(1))=zeros(size(W,1),1);
                density_list(i,1)=sum(sum(W,2)~=0);
                density_list(i,2)=sum(sum(W,1)~=0);
                density_list(i,3)=2;
            end
            density_list(i,4)=IndS(1);
            density_list(i,5)=dS;
            density_list(i,6)=IndT(1);
            density_list(i,7)=dT;
            %original density of W
            %W2=W.*25;
 
            %density_list(i,8)=(sum(sum(W,2)))/(sqrt(density_list(i,1)*density_list(i,2)))^lambda_w;
            density_list(i,8)=length(find(W>threshold))/(sqrt(density_list(i,1)*density_list(i,2)))^lambda_w;
        
            if density_list(i,2)>1
            
                %%% Construct u based on extracted edges 
                %%% first, collect removed edges up to current iteration i 
                remove_edges=[];
                for j=1:i
                    if density_list(j,3)==2
                        remove_edges =[remove_edges density_list(j,6)];
                    end
                end
    
                %%% now, collect retained edges 
                retain_edges0=setdiff(1:size(W,2),remove_edges); 
                retain_edges=t_in0(retain_edges0);
                edge_all=zeros(1,num_region*(num_region-1)/2);
                edge_all(retain_edges)=1;
                
                %%% construct correpsonding network u
                u=squareform(edge_all);
%                 retain_nodei=edge_node_list(retain_edges,2);
%                 retain_nodej=edge_node_list(retain_edges,3);
    
                %%% contruct u with entries= 0,1 
%                 u=zeros(num_region,num_region);  %e.g., 50*50
%                 for k=1:size(retain_edges,2)
%                     u(retain_nodei(k),retain_nodej(k))=1;
%                 end 
%                 u=u+u.';
                % figure;imagesc(u);ax=gca;ax.FontSize=18;
    
                [inlist outlist,max_density]=greedy_innerlayer(u,lambda_u);
                density_list(i,9)=max_density;  %density of U
%                 u_return = u([inlist outlist],[inlist outlist]);
%                 figure;imagesc(u_return);colormap jet;colorbar;
    
                %map detected voxels back to edges 
                edge_idx=[];
                for m=1:size(inlist,2)
                    for n=1:size(inlist,2)
                        if n>m
                            idx= find(edge_node_list(:,2)==inlist(m) & edge_node_list(:,3)==inlist(n));
                            edge_idx=[edge_idx idx];
                        end
                    end 
                end 
                 
                %find overlap between edge_idx and retain_edges
                up_retain_edges0=intersect(edge_idx, retain_edges); %update retained edges in original index
                
                %change to the indexing system based on shrinked W
                up_retain_edges=[];
                for p=1:length(up_retain_edges0)
                    up_retain_edges=[up_retain_edges find(t_in0==up_retain_edges0(p))];  
                end 
    
                retain_edge_list{i,1}=up_retain_edges;
    
                %% Update W based on u1(densly sub-graph extracted from u)
                %%% % by truncating origiinal edges
                W_trct=W(:,up_retain_edges); %note that rows are zeroed already, so selecting all rows are equivlent to selecting rows in current sub-network 
                %density_list(i,10)=(sum(sum(W_trct,2)))/(sqrt(density_list(i,1)*length(up_retain_edges)))^lambda_w;
                density_list(i,10)=length(find(W_trct>threshold))/(sqrt(density_list(i,1)*length(up_retain_edges)))^lambda_w;
            
                % add in the density of u   
                density_list(i,11)=density_list(i,9)+density_list(i,10);
            else 
                break
            end 
        
        end
    
    
        %ouput W
        [dST, indST]=max(density_list(:,11));
        indST
        dW=density_list(indST,10); %stopping criterion- density of truncated W

        %%%% SNPs:
        remove_s=[]; %SNPs that were removed    
        for q=1:indST
            if density_list(q,3)==1
                remove_s =[remove_s density_list(q,4)];
            end
        end
        retain_s=setdiff(1:size(W,1),remove_s);
        retain_s_vec{l}=retain_s;


        %%%% edges:
        retain_t=retain_edge_list{indST,1};
        retain_t_vec{l}=retain_t;

        dST_vec(l,1)=dST;
        dST_vec(l,2)=dW;
    end

    denWU_max=max(dST_vec(:,1));
    c_max = find(dST_vec == denWU_max);
    denW_max=dST_vec(c_max,2);
    s_in = retain_s_vec{c_max(1)};
    t_in = retain_t_vec{c_max(1)};
end
