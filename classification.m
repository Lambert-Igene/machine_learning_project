addpath("C:\Users\igene\OneDrive - clarkson.edu\emg\Features")
load features0.mat;
load features1.mat;
%load('T2data.mat');
%N1(:,end+1) = 
T=[N1;N];
%L=Data(:,end);
%Data=[L0;L1];
%L=Data(:,end);
%N
%[train,test] = crossvalind('HoldOut',L,.5);
%Data_train=Data(train,end-1);
%L_train=L(train,1);
%Data_test=Data(test,end-1);
%L_test=L(test,1);

%Md1 = fitcsvm(Data_train,L_train,'Standardize',true,'KernelFunction','RBF','KernelScale','auto'); % Bayes' Optimization ??.
%imp = predictorImportance(Mdl);
% test_accuracy = sum((predict(Md1,Data_test) == L_test))/length(L_test)*100;
% train_accuracy = sum((predict(Md1,Data_train) == L_train))/length(L_train)*100;



