
%% Load Data
load('../data/DataPartition3.mat')
Xbow_train = bagOfWords(gene_names_all(keep_g==1),X_train');
Xbow_test1 = bagOfWords(gene_names_all(keep_g==1),X_test1');
% Xtfidf = tfidf(bagOfWords(gene_names_all(keep_g==1),X'))
%% LDA
% length(unique(DataSet.celldata.cell_IDs)) Potential topics
% 
% numTopics = 73;
% TopicConcentration; %alpha
% % WordConcentration; %beta
% LDAmod.properties.Verbose = 1;
% LDAmod.properties.Solver = 'savb';
% LDAmod.properties.LogLikelihoodTolerance = 0.001;
% LDAmod.properties.FitTopicProbabilities = false;
% LDAmod.properties.FitTopicConcentration = false;

numTopicsRange = [2:2:20];

for i=1:4 %length(numTopicsRange)
    tic
    numTopics = numTopicsRange(i);
    mdlLDA{i} = fitlda(Xbow_train,numTopics,...
        'Verbose',1,...                             %display outputs
        'Solver','savb',...                         %solver used to optimize
        'LogLikelihoodTolerance',0.001,...
        'FitTopicProbabilities',true,...
        'FitTopicConcentration',true,...
        'ValidationData',Xbow_test1);
    toc
    [~,validationPerplexityLDA(i)] = logp(mdlLDA{i},Xbow_test1);
    toc
    timeElapsedLDA(i) = mdlLDA{i}.FitInfo.History.TimeSinceStart(end);
end

% learn rate decay 0.5
% 

figure; subplot(2,1,1)
bar(numTopicsRange,validationPerplexityLDA);
xlabel('Number of Topics');
ylabel('Model Preplexity');
title("LDA Preplexity")
subplot(2,1,2)
bar(numTopicsRange,timeElapsedLDA);
xlabel('Number of Topics');
ylabel('Training Time (s)');
title("LDA Training Time")

save('LDAmod_2.mat','mdlLDA')

%%

for i=1:length(mdlLDA)
%     LDAphi(i) = mdlLDA{i}.CorpusTopicProbabilities;
    LDAalpha(i) = mdlLDA{i}.TopicConcentration;
    LDAbeta(i) = mdlLDA{i}.WordConcentration;
end

figure
bar(numTopicsRange,LDAalpha)
xlabel('Number of Topics')
ylabel("Topic Concentration (alpha)")

%% LSA
% length(unique(DataSet.celldata.cell_IDs)) Potential topics
for i=1:length(numTopicsRange)
    tic
    numTopics = numTopicsRange(i);
    mdlLSA{i} = fitlsa(Xbow_train,numTopics,...
    'FeatureStrengthExponent',2);                           %default: 2
    timeElapsedLSA(i) = toc
end

% figure
% bar(numTopicsRange,validationPerplexityLSA);
% xlabel('Number of Topics');
% ylabel('Model Preplexity');

figure
bar(numTopicsRange,timeElapsedLSA);
xlabel('Number of Topics');
ylabel('Training Time');
title("LSA Training Time")


save('LSAmod_2.mat','mdlLSA')



%%
for i=1:length(numTopicsRange)
    ComponentWeights(i) = sum(mdlLSA{i}.ComponentWeights);
end

figure; subplot(2,1,1)
bar(ComponentWeights)
xlabel('Number of Topics');
ylabel('Component Weights');
title("LSA Component Sum")
subplot(2,1,2)
bar(numTopicsRange,timeElapsedLSA);
xlabel('Number of Topics');
ylabel('Training Time');
title("LSA Training Time")

%% NNMF
% 'Options',statset('UseParallel',true)
% parpool

% numTopicsRange = [10,20,40,73,4,100];
for i=1:length(numTopicsRange)
    tic
    numTopics = numTopicsRange(i);
    [mdlNNMF{i}.W,mdlNNMF{i}.H,mdlNNMF{i}.D] = nnmf(X_train,numTopics,...
        'Algorithm','als',...
        'Options',statset('UseParallel',true));
    Dnorm(i) = mdlNNMF{i}.D;
    mdlNNMF{i}.timeElapsedNNMF = toc
    timeElapsedNNMF(i) = mdlNNMF{i}.timeElapsedNNMF;
end

figure
bar(numTopicsRange,Dnorm);
xlabel('Number of Topics');
ylabel('RMS Residual');
title("NNMF Reconstruction Error")
figure
bar(numTopicsRange,timeElapsedNNMF);
xlabel('Number of Topics');
ylabel('Computation Time (s)');
title("NNMF Computation Time")


save('NMFmod_2.mat','mdlNNMF')


%% Export Topic Loadings

k = 30; %Top 30 genes in each topic
numTopicsRange = [2:2:20];

%LDA
for i=1:length(numTopicsRange)
    clear tbl
    for j=1:numTopicsRange(i)
        tbl_col = topkwords(mdlLDA{i},k,j);
        tbl(:,j) = [strcat('Topic',num2str(j)); tbl_col.Word];
    end
    filename = strcat("../results/loadings/LDA",num2str(numTopicsRange(i)),".txt");
    writematrix(tbl,filename,'Delimiter','tab');
%     writematrix(filename,"../results/filedirectory_load.txt",'WriteMode','append');
end

%LSA
for i=1:length(numTopicsRange)
    clear tbl
    for j=1:numTopicsRange(i)
        [~,Idx]=sort(mdlLSA{i}.WordScores(:,j),"descend");
        tbl(:,j) = [strcat('Topic',num2str(j)); string(gene_names_filtered(Idx(1:k)))];
    end
    filename = strcat("../results/loadings/LSA",num2str(numTopicsRange(i)),".txt");
    writematrix(tbl,filename,'Delimiter','tab');
%     writematrix(filename,"../results/filedirectory_load.txt",'WriteMode','append');
    
end

%NMF
for i=1:length(numTopicsRange)
    clear tbl
    for j=1:numTopicsRange(i)
        [~,Idx]=sort(mdlNNMF{i}.W(:,j),"descend");
        tbl(:,j) = [strcat('Topic',num2str(j)); string(gene_names_filtered(Idx(1:k)))];
    end
    filename = strcat("../results/loadings/NMF",num2str(numTopicsRange(i)),".txt");
    writematrix(tbl,filename,'Delimiter','tab');
%     writematrix(filename,"../results/filedirectory_load.txt",'WriteMode','append');
    
end

%% Export Document Loadings

%LDA
for i=1:length(numTopicsRange)
    filename = strcat("../results/matrices/LDA",num2str(numTopicsRange(i)),".txt");
    LDAmat = transform(mdlLDA{i},X');
    writematrix(LDAmat,filename,'Delimiter','tab')
    writematrix(filename,"../results/filedirectory_mtx.txt",'WriteMode','append');
end

%LSA
for i=1:length(numTopicsRange)
    filename = strcat("../results/matrices/LSA",num2str(numTopicsRange(i)),".txt");
    LSAmat = transform(mdlLSA{i},X');
    writematrix(LSAmat,filename,'Delimiter','tab');
    writematrix(filename,"../results/filedirectory_mtx.txt",'WriteMode','append');
end

%NNMF
for i=1:length(numTopicsRange)
    filename = strcat("../results/matrices/NMF",num2str(numTopicsRange(i)),".txt");
    NMFmat = X'*mdlNNMF{i}.W;
    writematrix(NMFmat,filename,'Delimiter','tab');
    writematrix(filename,"../results/filedirectory_mtx.txt",'WriteMode','append');
end

filename = strcat("../results/matrices/cellnames.txt");
writecell(CellLabels,filename,'Delimiter','tab');