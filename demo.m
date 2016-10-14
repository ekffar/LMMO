train_positive_set = 'training_positive.fasta';
train_negative_set = 'training_negative.fasta';
test_positive_set = 'test_positive.fasta';
test_negative_set = 'test_negative.fasta';

load Initial_Motif

[ Train_AUC, Test_AUC, PWM ] = LMMO_TCBB( train_positive_set, train_negative_set, test_positive_set, test_negative_set, Initial_Motif, 0.05);