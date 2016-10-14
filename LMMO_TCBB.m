function [ Train_AUC, Test_AUC, PWM ] = LMMO_TCBB( train_positive_set, train_negative_set, test_positive_set, test_negative_set, w0, p)
% LMMO for refining motifs.
% Input Arguments:
% train_positive_set, train_negative_set, test_positive_set,
% test_negative_set: character vectors that specify the file names (or
% paths and file names) of the positive/negative sequence sets for
% training/testing, respectively. ( Note that test_positive_set and
% test_negative_set are only used for evaluating the refined motifs,)
% The referenced files should be FASTA-formatted.
% p: the ratio of initially selected positive-negative sequence pairs for training.
% w0: the inpput PWM that needs to be refined.

% Output Arguments:
% Train_AUC: a vector which records the training AUCs of LMMO during each iteration.
% Test_AUC: a vector which records the training AUCs of LMMO during each iteration.
% PWM: the refined motif.

% Currently, this code is written such that it can run properly using the built-in functions of Matlab only. However, It is strongly recommended to also install Mosek before running LMMO,
% otherwise the LP in LMMO would require much longer time to compute. The details about how to install Mosek can be found at: http://docs.mosek.com/7.1/quickstart/Using_MOSEK_from_MATLAB.html 

motif_length = length(w0)/4;

[ P_Train_Overall_Data_Numerical_Sub, P_Tr_Bag_Size ] = orthogonal_coding( train_positive_set, motif_length );
[ N_Train_Overall_Data_Numerical_Sub, N_Tr_Bag_Size ] = orthogonal_coding( train_negative_set, motif_length );
[ P_Test_Overall_Data_Numerical_Sub, P_Te_Bag_Size ] = orthogonal_coding(test_positive_set, motif_length);
[ N_Test_Overall_Data_Numerical_Sub, N_Te_Bag_Size ] = orthogonal_coding(test_negative_set, motif_length);

Instance = [P_Train_Overall_Data_Numerical_Sub;N_Train_Overall_Data_Numerical_Sub];
Bag_Label = [ones(length(P_Tr_Bag_Size),1); -ones(length(N_Tr_Bag_Size),1)];

Test_Overall_Data_Numerical_Sub = [P_Test_Overall_Data_Numerical_Sub; N_Test_Overall_Data_Numerical_Sub];
Test_Data_Gnd = [ones(length(P_Te_Bag_Size),1); -ones(length(N_Te_Bag_Size),1)];
Tr_Bag_Size  = [P_Tr_Bag_Size; N_Tr_Bag_Size];
Te_Bag_Size = [P_Te_Bag_Size; N_Te_Bag_Size];

clear P_Train_Overall_Data_Numerical_Sub P_Tr_Bag_Size N_Train_Overall_Data_Numerical_Sub N_Tr_Bag_Size
clear P_Test_Overall_Data_Numerical_Sub P_Te_Bag_Size N_Test_Overall_Data_Numerical_Sub N_Te_Bag_Size

w0 = w0(:);

Initial_NN = 1; 
% Initial_NN: The initial number of top-scored subsequences that are
% seleced for each negative sequences.  In the experiments of the LMMO
% paper, this parameter is set as 1. 

dim = size(Instance,2); % The dimension of coding feature for each subsequence.
np = length(find(Bag_Label==1)); % The number of positive sequences.
nn = length(find(Bag_Label==-1)); % The number of negative sequences.

iter_num = 200;
%The maximum number of interations. It is deliberately set as a large number to ensure that LMMO does not terminate until the AUC value does not improve further:in practice, LMMO generally
%converged after less than 40 interations.

Select_Num = ceil(np*nn*p)+iter_num;

Tr_Instance_Bag_Index = cell(length(Bag_Label),1);
index_right = 0;
for i = 1:length(Tr_Bag_Size)
    index_left = index_right;
    index_right = index_right+Tr_Bag_Size(i);
    Tr_Instance_Bag_Index{i} = [index_left+1:index_right];
end

Te_Instance_Bag_Index = cell(length(Te_Bag_Size),1);
index_right = 0;
for i = 1:length(Te_Bag_Size)
    index_left = index_right;
    index_right = index_right+Te_Bag_Size(i);
    Te_Instance_Bag_Index{i} = [index_left+1:index_right];
end

Tmp_Tr_Instance_Bag_Index = Tr_Instance_Bag_Index;% During the iterations, Tmp_Tr_Instance_Bag_Index will record the index of subsequences that have not been added to the constraints of the LP.
Score = Instance*w0;

P_Index = zeros(np,1);
N_Index = zeros(nn*Initial_NN,1);

iA = 0;
iAeq = 0;

for i = 1:length(Bag_Label)
    index = Tr_Instance_Bag_Index{i};   
    
    if Bag_Label(i) == 1        
        [~, index_i] = max(Score(index));
        P_Index(iAeq+1) = index(index_i);
        iAeq = iAeq + 1;
    else
        [~,order] = sort(-Score(index));
        N_Index(iA+1:iA+Initial_NN,:) = index(order(1:Initial_NN));
        iA = iA + Initial_NN;   
        Tmp_Tr_Instance_Bag_Index{i}(order(1:Initial_NN)) = [];
    end  
end

Aeq = [Instance(P_Index,:) -speye(np) sparse(np,nn) sparse(np,Select_Num)];
beq = zeros(np,1);

A2 = zeros(nn*Initial_NN, nn);

z = 0;
c = 0;
for i = 1:np+nn  
    if Bag_Label(i) == -1
        c = c + 1;
        A2(z+1:z+Initial_NN,c) = -1;
        z = z + Initial_NN;
    end
end

A0 = Instance(N_Index,:);
A = [A0 sparse(nn*Initial_NN,np) sparse(A2) sparse(nn*Initial_NN,Select_Num)]; 
b = sparse(size(A,1),1);

Coef1 = [1:np];
Coef1 = Coef1';
Coef1 = repmat(Coef1,1,nn);
Coef1 = reshape(Coef1,1,[]);

Coef2 = np + [1:nn];
Coef2 = repmat(Coef2,np,1);
Coef2 = reshape(Coef2,1,[]);

ColID = [Coef1' Coef2'];

RowID = [1:Select_Num]';
RowID = [RowID RowID];

order = randperm(np*nn);
order = order(1:Select_Num);

A1 = sparse(reshape(RowID,[],1),reshape(ColID(order,:),[],1),reshape(repmat([-1 1],Select_Num,1),[],1),Select_Num,np+nn);

A1 = [sparse(Select_Num,dim) A1 -speye(Select_Num)];
b1 = -ones(Select_Num,1);

A0 = A;
A = [A; A1];
b = [b; b1];

Incre_A_Part_2 = [sparse(nn,np) -speye(nn) sparse(nn,Select_Num)];
Incre_A_Part_1 = zeros(nn,dim);

clear b1 A1  A2

Aeq = sparse(Aeq);
A = sparse(A);

% options = mskoptimset('Display','off','Diagnostics','off');

c_w = [];
for j = 1:iter_num
      
    w = w0;
    
    options.Display = 'off';
    
    disp(['Iteration: ' num2str(j)]);
    
    if j > 1   
        
        Score = Instance*w;
        z = 0;
        Bag_Score = zeros(np+nn,1);
        change_index = zeros(np,1);
        for i = 1:length(Bag_Label)
            index = Tmp_Tr_Instance_Bag_Index{i};
            [Bag_Score(i), index_i] = max(Score(index));
            if Bag_Label(i) == 1
                z = z + 1;
                change_index(z) = index(index_i);
            end
        end % Identify the optimal matching position for each positive sequence.
        
%         order = randperm(np*nn);  % Optional randomlization.
%         order = order(1:Select_Num); % Optional randomlization.
        
        RowID = [1:Select_Num - iter_num + j]';
        RowID = [RowID RowID];
        
        Tmp_ColID = ColID(order(1:Select_Num - iter_num + j - 1), :);
        
        [~, max_p_index] = min(Bag_Score(Bag_Label>0));
        [~, max_n_index] = max(Bag_Score(Bag_Label<0));
        max_n_index = max_n_index + np;
        
        Tmp_ColID = [Tmp_ColID; max_p_index max_n_index];
        
        A1 = sparse(reshape(RowID,[],1),reshape(Tmp_ColID,[],1),reshape(repmat([-1 1],size(RowID,1),1),[],1),Select_Num,np+nn);
        A1 = [A0; sparse(Select_Num,dim) A1 -speye(Select_Num)];
        A(1:size(A1,1),:) = [];
        A = [A1; A];
        
        a = Instance(change_index,:)*randn(dim,1);
        [~, ib, ic] = unique(a);
        mapping_ma = sparse(ic,[1:np]',ones(np,1),length(ib),np);
        Aeq = [sparse(Instance(change_index(ib),:)) -speye(length(ib)) sparse(length(ib),nn) sparse(length(ib),Select_Num)];
        beq = zeros(length(ib),1);        
        A_New = A(:,dim+1:dim+np);
        A_New = A_New*mapping_ma';
        A_New = [A(:,1:dim) A_New A(:,dim+np+1:end)];
        A_New = sparse(A_New); %Simplify the structure of the LP.
        
        f = [sparse(dim+length(ib)+nn,1);ones(Select_Num - iter_num + j,1);sparse(iter_num-j,1)];
        lb = [-Inf*ones(dim+length(ib)+nn,1);zeros(Select_Num,1)];
        lb = sparse(lb);
        tic
        c_w = linprog(f,A_New,b,Aeq,beq,lb,[],c_w,options);
        toc
        w = c_w(1:length(w0));
    end
    
    Score = Instance*w;
    
    Bag_Score = Bag_Label;
    
    z2 = 0;
    n_indicator = zeros(nn,1);

    for i = 1:np+nn
        if Bag_Label(i) == 1
        else
            index = Tmp_Tr_Instance_Bag_Index{i};
            [Bag_Score(i), index_i] = max(Score(index));
            z2 = z2 + 1;
            if Bag_Score(i) >=  max(Score(Tr_Instance_Bag_Index{i}))
                Incre_A_Part_1(z2,:) = Instance(index(index_i),:);
                n_indicator(z2) = 1; % n_indicator indicates whether the z2-th negative sequence constains a subsequence that violates the contraints
                Tmp_Tr_Instance_Bag_Index{i}(index_i) = [];
            end
        end
    end 
    
    sum(n_indicator)
    
    A = [A;Incre_A_Part_1(n_indicator>0,:) Incre_A_Part_2(n_indicator>0,:)];       
    b = [b; sparse(sum(n_indicator),1)]; % Update the inequality constraints.
    
    for i = 1:length(Bag_Label)
        index = Tr_Instance_Bag_Index{i};
        Bag_Score(i) = max(Score(index));
    end
    
    [auc] = roc_calculation( Bag_Score, Bag_Label );
    
    Train_AUC(j) = auc;
        
    w = reshape(w,[],1);    
    
    Test_Data_Score = Test_Data_Gnd;
    Value1 = Test_Overall_Data_Numerical_Sub*w;
    for i = 1:length(Test_Data_Gnd)
        Test_Data_Score(i) = max(Value1(Te_Instance_Bag_Index{i}));
    end
    
    [auc] = roc_calculation( Test_Data_Score, Test_Data_Gnd );
    Test_AUC(j) = auc;
    
    if j > 1
        if Train_AUC(j) > Train_AUC(j-1) - 0.005
            w0 = w;
        else
            Train_AUC(j) = Train_AUC(j-1);
            Test_AUC(j) = Test_AUC(j-1);
        end
    end
    
    [Train_AUC(j)  Test_AUC(j)]
    
    PWM = reshape(w,4,[]);;
    
end

