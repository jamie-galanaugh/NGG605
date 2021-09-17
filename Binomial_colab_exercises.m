%exercise 1
n=10;
prelease=0.2;
k=0:10;
probs = binopdf(k,n,prelease)
%exercise 2
n=14;
k=8;
prob1=binopdf(k,n,0.1)
prob2=binopdf(k,n,0.2)  
prob3=binopdf(k,n,0.3)
prob4=binopdf(k,n,0.4)  
prob5=binopdf(k,n,0.5)
prob6=binopdf(k,n,0.6)  
prob7=binopdf(k,n,0.7)
prob8=binopdf(k,n,0.8)  
prob9=binopdf(k,n,0.9)
prob10=binopdf(k,n,1)  
%prob of 0.1 = 0!!
%most likely 0.6 prob of release (prob=0.2066)

%exercise 3
n1=14;
n2=14;
k1=8;
k2=5;
prelease=0.1
prob1=binopdf(k1,n1,prelease)
prob2=binopdf(k2,n2,prelease)
totalprob=prob1*prob2
totallogprob=log(prob1)+log(prob2)
ps=(0:0.01:1.0)';  
nps=length(ps);
probs=binopdf(repmat([k1 k2], nps, 1),repmat([n1 n2], nps, 1),repmat(ps, 1, 2));
subplot(2,1,1); cla reset; hold on;
ylabel('Likelihood');
likelihoodFcn = prod(probs,2);
plot(ps,likelihoodFcn); 
maxLikelihood=max(likelihoodFcn);
plot(ps(likelihoodFcn==maxLikelihood), maxLikelihood, 'ko');
subplot(2,1,2);cla reset;hold on;
ylabel('Log-likelihood');
logLikelihoodFcn = sum(log(probs),2);
plot(ps, logLikelihoodFcn);
maxLogLikelihood = max(logLikelihoodFcn);
plot(ps(logLikelihoodFcn==maxLogLikelihood), maxLogLikelihood, 'ko');
n=14;
prelease=0.3;
ks=0:n;
ps=0:0.01:1.0;
pdat=zeros(length(ps), 2);
TINY=0.0001;
for sampleSize = round(logspace(0,3,30))
    sCounts = binornd(n,prelease,sampleSize,1);
    subplot(3,1,1); cla reset; hold on;
    title(sprintf('Sample size = %d', sampleSize))
    ylabel('Probability');
    xlabel('Release count');
    xlim([ks(1) ks(end)]);
    sCountHistogram = hist(sCounts, ks);
    bar(ks, sCountHistogram./sum(sCountHistogram));
    plot(ks, binopdf(ks, n, prelease), 'ro-');
     for pp = 1:length(ps)
         probs = binopdf(sCounts, n, ps(pp));
         probs(probs<TINY) = TINY;
         pdat(pp,:) = [prod(probs), sum(log(probs))];
     end
     [phat, pci] = binofit(sum(sCounts),sampleSize*n);
      subplot(3,1,2); cla reset; hold on;
      ylabel('likelihood');
       plot(ps, pdat(:,1));
       maxp = max(pdat(:,1));
       plot(prelease.*[1 1], [0 maxp], 'r--');
       plot(pci, maxp.*[1 1], 'm-', 'LineWidth', 2);
       plot(phat, maxp, 'm*');
       plot(ps(pdat(:,1)==maxp), maxp, 'ko', 'MarkerSize', 12); 
       subplot(3,1,3); cla reset; hold on;
       ylabel('log likelihood');
       xlabel('Release probability');
       axis([0 1 log(TINY)*1000 0]);
       plot(ps, pdat(:,2));
       maxp = max(pdat(:,2));
       plot(prelease.*[1 1], [0 min(pdat(:,2))], 'r--');
       plot(pci, maxp.*[1 1], 'm-', 'LineWidth', 2);
       plot(phat, maxp, 'm*');
       plot(ps(pdat(:,2)==maxp), maxp, 'g*');
       pause(0.1);
end
%exercise 4
counts = [0	0 3 10 19 26 16 16 5 5 0 0 0 0 0]; 
n = length(counts)-1;   
ks = 0:n;               
nks = length(ks);       
ps = (0:0.01:1.0)';     
nps = length(ps);       

probs = binopdf( ...
   repmat(ks, nps, 1), ...    
   n, ...                     
   repmat(ps, 1, nks));       
countsMatrix = repmat(counts, nps, 1);
likelihoodFcn = prod(probs.^countsMatrix,2);
pHat_fromLiklihood = ps(likelihoodFcn==max(likelihoodFcn))
logLikelihoodFcn = sum(log(probs).*countsMatrix,2);
pHat_fromLogLikelihood = ps(logLikelihoodFcn==max(logLikelihoodFcn))
pHat = sum(counts.*ks)/(sum(counts)*n)

%% Exercise 5
n = 14;                       
k = 7;                        
pHat = binofit(k,n)           
pNull = 0.3;                  

1 - binocdf(k-1, n, pNull) 

