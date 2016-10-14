function [auc] = roc_calculation( b, label )
% 
% size(b)
% size(label)

N = b(label == -1);
P = b(label == 1);

Diff = repmat(P,1,length(N)) - repmat(N',length(P),1);

Diff(Diff>0) = 1;

Diff(Diff==0) = 0.5;

Diff(Diff<0) = 0;

auc = sum(Diff(:))/size(Diff,1)/size(Diff,2);

% save Diff Diff



% std(sum(Diff,2))


% 
% a = quantile(b(:),[0.001:0.001:1]);
% 
% % a = quantile([P;N],[0.005:0.005:1]);
% 
% % a = [0 a max(b)+ 2*eps];
% 
% % a
% na = length(a);
% 
% TP = zeros(na,1);
% TN = zeros(na,1);
% FP = zeros(na,1);
% FN = zeros(na,1);
% 
% for i = 1:na
%     TP(i) = length(P(P<a(i)))/length(P);
%     TN(i) = length(N(N>=a(i)))/length(N);
%     FN(i) = 1 - TP(i);
%     FP(i) = 1 - TN(i); 
% end
% 
% x1 = 1 - TN./(TN + FP + eps );
% x2 = TP./(TP + FN + eps );
% auc = 0;
% for i = 2:na
%     auc = (x1(i) - x1(i-1))*(x2(i) + x2(i-1))/2 + auc;
% end


% x1
% x2