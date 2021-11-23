function [ouputs]=zcoreheatmap3(datafluo,filename,timealig,name);

dfofCorr_sub2=datafluo{1, 3};

startnormal=timealig{3, 5};% add an if to make sure eliminate all the nan for all trials 
endnormal=timealig{7, 5};
completionnormal=timealig{9, 5}

startimpul=timealig{3, 6};
endimpul=timealig{7, 6};
  

startcompul=timealig{3, 7};
endcompul=timealig{7, 7};
completioncomp=timealig{9, 7};

startG=timealig{11,5};
endG=timealig{15,5};
completionG=timealig{17,5};


startSP=timealig{19,5};
endSP=timealig{23,5};
completionSP=timealig{25,5};



times=(-5:0.025:25.05);

% startnormal(find(isnan(startnormal)))=[];
% endnormal(find(isnan(endnormal)))=[];
% completionnormal(find(isnan(completionnormal)))=[];
% 
% 
% startimpul(find(isnan(startimpul)))=[];
% endimpul(find(isnan(endimpul)))=[]; %to eliminate nan for the matrix
% 
% 
% startcompul(find(isnan(startcompul)))=[];
% endcompul(find(isnan(endcompul)))=[];
% completioncomp(find(isnan(completioncomp)))=[];
% 
% 
% startG(find(isnan(startG)))=[];
% endG(find(isnan(endG)))=[];
% completionG(find(isnan(completionG)))=[];
% 
% 
% startSP(find(isnan(startSP)))=[];
% endSP(find(isnan(endSP)))=[];
% completionSP(find(isnan(completionSP)))=[];


%% graph start complete
if length(startnormal)>2
sampleci=(startnormal*40); samplec2i=round(sampleci);samplec2masi=samplec2i+200;samplec2menosi=samplec2i-200;

for k=1:(length(startnormal)-1)
    
    if startnormal(k,1)>200 
 trialscompl(k,:)= dfofCorr_sub2 (samplec2menosi(k,1):samplec2masi(k,1),1);
 zsctrialscomp(k,:)=zscore(trialscompl(k,:)); 
    end
end

% end completed first
sampleover=(endnormal *40); samplec2over=round(sampleover);samplec2masover=sampleover+200;samplec2mensoover=sampleover-200;
%dfofCorr_sub2i=dfofCorr_sub';

for k=1:(length(endnormal)-1)
    if endnormal(k,1)>200
 trialsendompled(k,:)= dfofCorr_sub2 (samplec2mensoover(k,1):samplec2masover(k,1),1);
 zsctrialendcompled(k,:)=zscore(trialsendompled(k,:)); 
  newcomp1(k,:)=[trialscompl(k,:) trialsendompled(k,:)];
  zscnewcomp1(k,:)=zscore(newcomp1(k,:));
    end
end

% end completed first
samplecomp=(completionnormal *40); samplec2comp=round(samplecomp);samplec2mascomp=samplecomp+200;samplec2menscomp=samplecomp-200;

%
for k=1:(length(completionnormal)-1)
    if completionnormal(k,1)>200
 trialsendompled1(k,:)= dfofCorr_sub2 (samplec2menscomp(k,1):samplec2mascomp(k,1),1);
 zsctrialendcompled(k,:)=zscore(trialsendompled(k,:)); 
  newcomp(k,:)=[trialscompl(k,:) trialsendompled(k,:) trialsendompled1(k,:)];
  zscnewcomp(k,:)=zscore(newcomp(k,:));
  
    end
end

bbbb=find(zscnewcomp(:,1)~=0);
tempmean=(sum(zscnewcomp))/(length(bbbb));%tempmean=smooth(tempmean,10,'moving');tempmean=tempmean';


figure(1)

subplot(2,3,1)
title('completed start')
image ( zscnewcomp,'CDataMapping','scaled');
xline(200,'-',{'first press'});
ylabel('trials');
colormap(jet)
axis ([0 400 1 inf]) 
caxis([-2 2]);
subplot(2,3,4)
plot(times,tempmean)
hold on
axis ([-5 5 -2 2]) 
xline(0,'-',{'first press'});
ylabel('zcore');xlabel('time sec')
subplot(2,3,2)
title('completed end')
image ( zscnewcomp,'CDataMapping','scaled');
ylabel('trials');
axis ([400 800 1 inf]) 
caxis([-2 2]);
xline(600,'-',{'last press'});
colormap(jet)
subplot(2,3,5)
plot(times,tempmean);
axis ([5 15 -2 2]) 
xline(10,'-',{'last press'});
 ylabel('zcore');xlabel('time sec')

subplot(2,3,3)
title('completed completion')
image ( zscnewcomp,'CDataMapping','scaled');
axis ([800 1200 1 inf]) 
xline(1000,'-',{'completion lever'});
caxis([-2 2]);
colormap(jet)
subplot(2,3,6)
plot(times,tempmean);
 xline(20,'-',{'completion lever'});
axis ([15 25 -2 2]) 
 ylabel('zcore');xlabel('time sec')


clear sampleci samplec2i samplec2masi samplec2menosi sampleover samplec2over samplec2masover samplec2mensoover



%% if you modify the time before or after you should modify this numbers for every loop

zscnewcomp2=tempmean;
aameanStart=trapz(zscnewcomp2(1,120:200));
aameanexec=trapz(zscnewcomp2(1,200:280));
aameanexec2=trapz(zscnewcomp2(1,520:600));
aameanend=trapz(zscnewcomp2(1,600:680));
aacomplet=trapz(zscnewcomp2(1,720:800));
aacomplet2=trapz(zscnewcomp2(1,800:880));


aameanStart1=trapz(zscnewcomp2(1,160:200));
aameanexec1=trapz(zscnewcomp2(1,200:240));
aameanexec21=trapz(zscnewcomp2(1,560:600));
aameanend1=trapz(zscnewcomp2(1,600:640));
aacomplet1=trapz(zscnewcomp2(1,760:800));
aacomplet21=trapz(zscnewcomp2(1,800:840));

ouputs(1,1)=aameanStart;
ouputs(1,2)=aameanexec;
ouputs(1,3)=aameanexec2;
ouputs(1,4)=aameanend;
ouputs(1,5)=aacomplet;
ouputs(1,6)=aacomplet2;

ouputs(9,1)=aameanStart1;
ouputs(9,2)=aameanexec1;
ouputs(9,3)=aameanexec21;
ouputs(9,4)=aameanend1;
ouputs(9,5)=aacomplet1;
ouputs(9,6)=aacomplet21;

[Where] = Where_file(filename);
SaveName= [strcat(Where,'figure_completed_fluorece_trials',name)];
     savefig(SaveName);
end
%%  graph imp

if length(startimpul)>2
sampleci=(startimpul*40); samplec2i=round(sampleci);samplec2masi=samplec2i+200;samplec2menosi=samplec2i-200;

for k=1:(length(startimpul)-1)
    
    if startimpul(k,1)>200
 trialsimp(k,:)= dfofCorr_sub2 (samplec2menosi(k,1):samplec2masi(k,1),1);
 zsctrialsimpuli(k,:)=zscore(trialsimp(k,:)); 
    end
end

% end completed first
sampleover=(endimpul *40); samplec2over=round(sampleover);samplec2masover=sampleover+200;samplec2mensoover=sampleover-200;
%dfofCorr_sub2i=dfofCorr_sub';

for k=1:(length(endimpul)-1)
    if endimpul(k,1)>200
 trialsendimp(k,:)= dfofCorr_sub2 (samplec2mensoover(k,1):samplec2masover(k,1),1);
 zsctrialendimp(k,:)=zscore(trialsendimp(k,:)); 
  newimp(k,:)=[trialsimp(k,:) trialsendimp(k,:)];
  zscnewimp(k,:)=zscore(newimp(k,:));
    end
end
times2=(-5:0.025:15.025);
bbbb2=find(zscnewimp(:,1)~=0);
tempmean2=(sum(zscnewimp))/(length(bbbb2));

figure(2)

subplot(2,3,1)
title('impulsive')
image ( zscnewimp,'CDataMapping','scaled');
colormap(jet)
axis ([0 400 1 inf]) 
xline(200,'-',{'first press'});
caxis([-2 2]);
subplot(2,3,4)
plot(times2,tempmean2)
axis ([-5 5 -2 2]) 
xline(0,'-',{'first press'});
ylabel('zcore');xlabel('time sec')
subplot(2,3,2)
image ( zscnewimp,'CDataMapping','scaled');
axis ([400 800 1 inf]);caxis([-2 2]);
xline(600,'-',{'last press'});

colormap(jet)
subplot(2,3,5)
plot(times2,tempmean2);
ylabel('zcore');xlabel('time sec')
axis ([5 15 -2 2])
xline(10,'-',{'last press'});

%% if you modify the time before or after you should modify this numbers for every loop


zscnewimp2=tempmean2;
bbmeanStart2=trapz(zscnewimp2(1,120:200));
bbmeanexec=trapz(zscnewimp2(1,200:280));
bbmeanexec2=trapz(zscnewimp2(1,520:600));
bbmeanend2=trapz(zscnewimp2(1,600:680));
% bbcomplet=trapz(zscnewimp2(1,720:800));
% bbcomplet2=trapz(zscnewimp2(1,800:880));

bbmeanStart21=trapz(zscnewimp2(1,160:200));
bbmeanexec1=trapz(zscnewimp2(1,200:240));
bbmeanexec21=trapz(zscnewimp2(1,560:600));
bbmeanend21=trapz(zscnewimp2(1,600:640));

clear sampleci samplec2i samplec2masi samplec2menosi sampleover samplec2over samplec2masover samplec2mensoover

[Where] = Where_file(filename);
SaveName= [strcat(Where,'figure__fluore_impulsive_trials',name)];
     savefig(SaveName);
ouputs(2,1)=bbmeanStart2;
ouputs(2,2)=bbmeanexec;
ouputs(2,3)=bbmeanexec2;
ouputs(2,4)=bbmeanend2;
% outputs(2,5)=bbcomplet;
% outputs(3,6)=bbcomplet2;

ouputs(10,1)=bbmeanStart21;
ouputs(10,2)=bbmeanexec1;
ouputs(10,3)=bbmeanexec21;
ouputs(10,4)=bbmeanend21;

end

%% compulsive

if length(startcompul)>2

sampleci=(startcompul*40); samplec2i=round(sampleci);samplec2masi=samplec2i+200;samplec2menosi=samplec2i-200;

for k=1:(length(startcompul)-1)
    
    if startcompul(k,1)>200
 trialscompul(k,:)= dfofCorr_sub2 (samplec2menosi(k,1):samplec2masi(k,1),1);
 zsctrialscompul(k,:)=zscore(trialscompul(k,:)); 
    end
end

% end completed first
sampleover=(endcompul *40); samplec2over=round(sampleover);samplec2masover=sampleover+200;samplec2mensoover=sampleover-200;
%dfofCorr_sub2i=dfofCorr_sub';

for k=1:(length(endcompul)-1)
    if endcompul(k,1)>200
 trialsendcompul(k,:)= dfofCorr_sub2 (samplec2mensoover(k,1):samplec2masover(k,1),1);
 zsctrialendcompul(k,:)=zscore(trialsendcompul(k,:)); 
  newcompu2(k,:)=[trialscompul(k,:) trialsendcompul(k,:)];
  zscnewcompu2(k,:)=zscore(newcompu2(k,:));
    end
end

% end completed first
samplecompcompul=(completioncomp *40); samplec2compcompul=round(samplecompcompul);samplec2mascompcompul=samplecompcompul+200;samplec2menscompcompul=samplecompcompul-200;

%
for k=1:(length(completioncomp)-1)
    if completioncomp(k,1)>200
 trialscompompled(k,:)= dfofCorr_sub2 (samplec2menscompcompul(k,1):samplec2mascompcompul(k,1),1);
 zsctrialendcompled(k,:)=zscore(trialscompompled(k,:)); 
  newcompul(k,:)=[trialscompul(k,:) trialsendcompul(k,:) trialscompompled(k,:)];
  zscnewcompul(k,:)=zscore(newcompul(k,:));
    end
end

bbbb3=find(zscnewcompul(:,1)~=0);
tempmean3=(sum(zscnewcompul))/(length(bbbb3));
figure(3)

subplot(2,3,1)
title('compulsive start')
image ( zscnewcompul,'CDataMapping','scaled');
colormap(jet)
axis ([0 400 1 inf]);caxis([-2 2]); 
xline(200,'-',{'first press'});

subplot(2,3,4)
plot(times,tempmean3)
axis ([-5 5 -2 2]) 
xline(0,'-',{'first press'});

ylabel('zcore');xlabel('time sec')

subplot(2,3,2)
title('compulsive end')
image ( zscnewcompul,'CDataMapping','scaled');
axis ([400 800 1 inf]);caxis([-2 2]); 
xline(600,'-',{'last press'});

colormap(jet)
subplot(2,3,5)
plot(times,tempmean3);
ylabel('zcore');xlabel('time sec')
xline(10,'-',{'first press'});

axis ([5 15 -2 2]) 


subplot(2,3,3)
title('compulsive completion')
image ( zscnewcompul,'CDataMapping','scaled');
axis ([800 1200 1 inf]) ;caxis([-2 2]);
colormap(jet)
xline(1000,'-',{'completion press'});

subplot(2,3,6)
plot(times,tempmean3);
axis ([15 25 -2 2])
ylabel('zcore');xlabel('time sec')
xline(20,'-',{'completion press'});
clear sampleci samplec2i samplec2masi samplec2menosi sampleover samplec2over samplec2masover samplec2mensoover

%% if you modify the time before or after you should modify this numbers for every loop
zscnewcompul2=tempmean3;
ccmeanStart2=trapz(zscnewcompul2(1,120:200));
ccmeanexec=trapz(zscnewcompul2(1,200:280));
ccmeanexec2=trapz(zscnewcompul2(1,520:600));
ccmeanend2=trapz(zscnewcompul2(1,600:680));
ccmeancompcomp=trapz(zscnewcompul2(1,720:800));
ccmeancompcomp2=trapz(zscnewcompul2(1,800:880));


ccmeanStart21=trapz(zscnewcompul2(1,160:200));
ccmeanexec1=trapz(zscnewcompul2(1,200:240));
ccmeanexec21=trapz(zscnewcompul2(1,560:600));
ccmeanend21=trapz(zscnewcompul2(1,600:640));
ccmeancompcomp1=trapz(zscnewcompul2(1,760:800));
ccmeancompcomp21=trapz(zscnewcompul2(1,800:840));


[Where] = Where_file(filename);
SaveName= [strcat(Where,'figure_fluore_compulsive_trials',name)];
     savefig(SaveName);

ouputs(3,1)=ccmeanStart2;
ouputs(3,2)=ccmeanexec;
ouputs(3,3)=ccmeanexec2;
ouputs(3,4)=ccmeanend2;
ouputs(3,5)=ccmeancompcomp;
ouputs(3,6)=ccmeancompcomp2;

ouputs(11,1)=ccmeanStart21;
ouputs(11,2)=ccmeanexec1;
ouputs(11,3)=ccmeanexec21;
ouputs(11,4)=ccmeanend21;
ouputs(11,5)=ccmeancompcomp1;
ouputs(11,6)=ccmeancompcomp21;

end 
%% graph guided complete
if length(startG)>2
sampleciG=(startG*40); samplec2iG=round(sampleciG);samplec2masiG=samplec2iG+200;samplec2menosiG=samplec2iG-200;

for k=1:(length(startG)-1)
    
    if startG(k,1)>200 
 trialscomplG(k,:)= dfofCorr_sub2 (samplec2menosiG(k,1):samplec2masiG(k,1),1);
 zsctrialscompG(k,:)=zscore(trialscomplG(k,:)); 
    end
end

% end completed first
sampleoverG=(endG *40); samplec2overG=round(sampleoverG);samplec2masoverG=sampleoverG+200;samplec2mensooverG=sampleoverG-200;
%dfofCorr_sub2i=dfofCorr_sub';

for k=1:(length(endG)-1)
    if endG(k,1)>200
 trialsendompledG(k,:)= dfofCorr_sub2 (samplec2mensooverG(k,1):samplec2masoverG(k,1),1);
 zsctrialendcompledG(k,:)=zscore(trialsendompledG(k,:)); 
  newcomp1G(k,:)=[trialscomplG(k,:) trialsendompledG(k,:)];
  zscnewcomp1G(k,:)=zscore(newcomp1(k,:));
    end
end

% end completed first
samplecompG=(completionG *40); samplec2compG=round(samplecompG);samplec2mascompG=samplecompG+200;samplec2menscompG=samplecompG-200;

%
for k=1:(length(endG)-1)
    if completionG(k,1)>200
 trialsendompled2G(k,:)= dfofCorr_sub2 (samplec2menscompG(k,1):samplec2mascompG(k,1),1);
 zsctrialendcompled(k,:)=zscore(trialsendompled2G(k,:)); 
  newcompG(k,:)=[trialscomplG(k,:) trialsendompledG(k,:) trialsendompled2G(k,:)];
  zscnewcompG(k,:)=zscore(newcompG(k,:));
    end
end

bbbb5=find(zscnewcompG(:,1)~=0);
tempmean5=(sum(zscnewcompG))/(length(bbbb5));
figure(5)

subplot(2,3,1)
title('completed start')
image ( zscnewcompG,'CDataMapping','scaled');
colormap(jet)
axis ([0 400 1 inf]);caxis([-2 2]); 
xline(200,'-',{'first press'});

subplot(2,3,4)
plot(times,tempmean5)
axis ([-5 5 -2 2]) 
ylabel('zcore');xlabel('time sec')
xline(0,'-',{'first press'});

subplot(2,3,2)
title('completed end')
image ( zscnewcompG,'CDataMapping','scaled');
axis ([400 800 1 inf]) ;caxis([-2 2]);
xline(600,'-',{'end press'});
colormap(jet)
subplot(2,3,5)
plot(times,tempmean5);
axis ([5 15 -2 2]) 
xline(10,'-',{'end press'});

ylabel('zcore');xlabel('time sec')
subplot(2,3,3)
title('completed completionG')
image ( zscnewcompG,'CDataMapping','scaled');
axis ([800 1200 1 inf]) ;caxis([-2 2]);
xline(1000,'-',{'completion press'});
colormap(jet)
subplot(2,3,6)
plot(times,tempmean5);
xline(20,'-',{'completion press'});

axis ([15 25 -2 2]) 


clear sampleciG samplec2iG samplec2masiG samplec2menosiG sampleoverG samplec2overG samplec2masoverG samplec2mensooverG

zscnewcomplG=tempmean5;
ddmeanStart2=trapz(zscnewcomplG(1,120:200));
ddmeanexec=trapz(zscnewcomplG(1,200:280));
ddmeanexec2=trapz(zscnewcomplG(1,520:600));
ddmeanend2=trapz(zscnewcomplG(1,600:680));
ddmeancompcomp=trapz(zscnewcomplG(1,720:800));
ddmeancompcomp2=trapz(zscnewcomplG(1,800:880));



ddmeanStart21=trapz(zscnewcomplG(1,160:200));
ddmeanexec1=trapz(zscnewcomplG(1,200:240));
ddmeanexec21=trapz(zscnewcomplG(1,560:600));
ddmeanend21=trapz(zscnewcomplG(1,600:640));
ddmeancompcomp1=trapz(zscnewcomplG(1,760:800));
ddmeancompcomp21=trapz(zscnewcomplG(1,800:840));

ouputs(4,1)=ddmeanStart2;
ouputs(4,2)=ddmeanexec;
ouputs(4,3)=ddmeanexec2;
ouputs(4,4)=ddmeanend2;
ouputs(4,5)=ddmeancompcomp;
ouputs(4,6)=ddmeancompcomp2;


ouputs(12,1)=ddmeanStart21;
ouputs(12,2)=ddmeanexec1;
ouputs(12,3)=ddmeanexec21;
ouputs(12,4)=ddmeanend21;
ouputs(12,5)=ddmeancompcomp1;
ouputs(12,6)=ddmeancompcomp21;

[Where] = Where_file(filename);
SaveName= [strcat(Where,'figure_fluore_comppleted_guided',name)];
     savefig(SaveName);

end
%% SP

if length(startSP)>2
sampleciSP=(startSP*40); samplec2iSP=round(sampleciSP);samplec2masiSP=samplec2iSP+200;samplec2menosiSP=samplec2iSP-200;

for k=1:(length(startSP)-1)
    
    if startSP(k,1)>200 
 trialscomplSP(k,:)= dfofCorr_sub2 (samplec2menosiSP(k,1):samplec2masiSP(k,1),1);
 zsctrialscompSP(k,:)=zscore(trialscomplSP(k,:)); 
    end
end

% end completed first
sampleoverSP=(endSP *40); samplec2overSP=round(sampleoverSP);samplec2masoverSP=sampleoverSP+200;samplec2mensooverSP=sampleoverSP-200;
%dfofCorr_sub2i=dfofCorr_sub';

for k=1:(length(endSP)-1)
    if endSP(k,1)>200
 trialsendompledSP(k,:)= dfofCorr_sub2 (samplec2mensooverSP(k,1):samplec2masoverSP(k,1),1);
 zsctrialendcompledSP(k,:)=zscore(trialsendompledSP(k,:)); 
  newcomp1SP(k,:)=[trialscomplSP(k,:) trialsendompledSP(k,:)];
  zscnewcomp1SP(k,:)=zscore(newcomp1SP(k,:));
    end
end

% end completed first
samplecompSP=(completionSP *40); samplec2compSP=round(samplecompSP);samplec2mascompSP=samplecompSP+200;samplec2menscompSP=samplecompSP-200;

%
for k=1:(length(endSP)-1)
    if completionSP(k,1)>200
 trialsendompled2SP(k,:)= dfofCorr_sub2 (samplec2menscompSP(k,1):samplec2mascompSP(k,1),1);
 zsctrialendcompled(k,:)=zscore(trialsendompled2SP(k,:)); 
  newcompSP(k,:)=[trialscomplSP(k,:) trialsendompledSP(k,:) trialsendompled2SP(k,:)];
  zscnewcompSP(k,:)=zscore(newcompSP(k,:));
    end
end

bbbb6=find(zscnewcompSP(:,1)~=0);
tempmean6=(sum(zscnewcompSP))/(length(bbbb6));


figure(6)

subplot(2,3,1)
title('completed start')
image ( zscnewcompSP,'CDataMapping','scaled');
colormap(jet)
axis ([0 400 1 inf]);caxis([-2 2]); 
xline(200,'-',{'first press'});
subplot(2,3,4)
plot(times,tempmean6)
axis ([-5 5 -2 2]) 
xline(0,'-',{'first press'});
ylabel('zcore');xlabel('time sec')
subplot(2,3,2)
title('completed end')
image ( zscnewcompSP,'CDataMapping','scaled');
axis ([400 800 1 inf]);caxis([-2 2]);
xline(600,'-',{'last press'});
colormap(jet)
subplot(2,3,5)
plot(times,tempmean6);
xline(10,'-',{'last press'});
axis ([5 15 -2 2]) 
subplot(2,3,3)
title('completed completionG')
image ( zscnewcompSP,'CDataMapping','scaled');
axis ([800 1200 1 inf]);caxis([-2 2]); 
colormap(jet)
xline(1000,'-',{'completion press'});

subplot(2,3,6)
plot(times,tempmean6);
xline(20,'-',{'completion press'});
axis ([15 25 -2 2]) 


clear sampleciSP samplec2iSP samplec2masiSP samplec2menosiSP sampleoverSP samplec2overSP samplec2masoverSP samplec2mensooverSP


zscnewcomplSP=tempmean6;
eemeanStart2=trapz(zscnewcomplSP(1,120:200));
eemeanexec=trapz(zscnewcomplSP(1,200:280));
eemeanexec2=trapz(zscnewcomplSP(1,520:600));
eemeanend2=trapz(zscnewcomplSP(1,600:680));
eemeancompcomp=trapz(zscnewcomplSP(1,720:800));
eemeancompcomp2=trapz(zscnewcomplSP(1,800:880));



eemeanStart21=trapz(zscnewcomplSP(1,160:200));
eemeanexec1=trapz(zscnewcomplSP(1,200:240));
eemeanexec21=trapz(zscnewcomplSP(1,560:600));
eemeanend21=trapz(zscnewcomplSP(1,600:640));
eemeancompcomp1=trapz(zscnewcomplSP(1,760:800));
eemeancompcomp21=trapz(zscnewcomplSP(1,800:840));





ouputs(5,1)=eemeanStart2;
ouputs(5,2)=eemeanexec;
ouputs(5,3)=eemeanexec2;
ouputs(5,4)=eemeanend2;
ouputs(5,5)=eemeancompcomp;
ouputs(5,6)=eemeancompcomp2;

ouputs(13,1)=eemeanStart21;
ouputs(13,2)=eemeanexec1;
ouputs(13,3)=eemeanexec21;
ouputs(13,4)=eemeanend21;
ouputs(13,5)=eemeancompcomp1;
ouputs(13,6)=eemeancompcomp21;

[Where] = Where_file(filename);
SaveName= [strcat(Where,'figure_fluore_comppleted_self_paced',name)];
     savefig(SaveName);

end

end
%%








