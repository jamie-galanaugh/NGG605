 function [datas3]=completion_Task_4(filename,name,timealinfotometry);

%%salidas=primera es datos de secuencia 1 y el 2 datoss de secuencia 2
%%incluye datos de todas las los tipos de cadenas
% el segundo renglon tiene 2 cosas sobre los guiados una son los parametros
%ver funcion y la otra son los indices de los ensayos que menciona
%se aplica en el tercer renglon lo mismo para los no guiados 
% PARAMETROS 1transicion 2tiempo latenciahead 3duracionS1 4duracionS2 5presionesS1 6presionesS2 7velocidadS1 8velocidadS2   

     %% datoss para sort
     %1 tiempo transicion
     %2 es para latencia conheaentry
     %3 tiempo duracion seq 1
     %4 tiempo duracion seq2
     %5 presiones S1
     %6 presiones S2
     %7 Velocidad P1
     %8 velocidad P2


%%arrastrar los datos de conducta y crear variable en conducta

[conductas2]=leertxt(filename);
%% colocar el tiempo de analisis la diferencia entre la primera pres y la del nev



%%separar la conducta y generar los tiempos
% sizes=size(conducta);
% volt=conducta';
% voltss=reshape(volt,1,5*length(conducta));
volts=conductas2;
volts(:,2)=fix(volts);%tiempo normal
volts(:,3)=volts(:,1)-volts(:,2);%banderas
volts(:,4)=(cumsum(volts(:,2)))/100;%tiempo acumulado seg

%%
%banderas
for i=1:length(volts);
    licks(i,1)=(volts(i,3)>0.109 && volts(i,3)<0.115);
    headentry(i,1)=(volts(i,3)>0.119 && volts(i,3)<0.121);
    press1(i,1)= (volts(i,3)>0.09 && volts(i,3)<0.105);
    press2(i,1)= (volts(i,3)>0.69 && volts(i,3)<0.705);
    reward(i,1)=(volts(i,3)>0.19 && volts(i,3)<0.21);
    rewardG(i,1)=(volts(i,3)>0.129 && volts(i,3)<0.131);
    rewardNG (i,1)=(volts(i,3)>0.139 && volts(i,3)<0.141);
    castigoP1G(i,1)=(volts(i,3)>0.549 && volts(i,3)<0.551);
   % castigoP2G(i,1)=(volts(i,3)>0.659 && volts(i,3)<0.661);
    castigoP1NG(i,1)=(volts(i,3)>0.579 && volts(i,3)<0.581);
   % castigoP2NG(i,1)=(volts(i,3)>0.679 && volts(i,3)<0.681);
end

tiempoP1parachecarnev=volts(press2,end);
salidas{6,1}=tiempoP1parachecarnev;

%% NUEVAS MARCAS PARA BANDERAS DEBIDO A QUE NO SON UN NUMERO CONSTANTE EN EL MEDPC
licks=licks*1;
headentry=headentry*1;
press1=press1*3;%ojo ver con que palanca empieza la diferencia de palancas tiene que er -1
press2=press2*4;%ojo ver con que palanca empieza 
reward=reward*5;
rewardG=rewardG*6;
rewardNG=rewardNG*7;
castigoP1G=castigoP1G*1;
%castigoP2G=castigoP2G*1;
castigoP1NG=castigoP1NG*1;
%castigoP2NG=castigoP2NG*1;
findpress2=find(press1);
findpress1=find(press2);

tiempo_acum=volts(:,end);tiempo_acum=tiempo_acum;
timespress2=tiempo_acum(findpress2,1);
timespressdiff=diff(timespress2);
timeleverP1bb=tiempo_acum(findpress1,1);
timespress1diff=diff(timeleverP1bb);


lickyhead=(headentry+castigoP1G+castigoP1NG);
banderas_mod=licks+headentry+press1+press2+reward+rewardG+rewardNG+castigoP1G+castigoP1NG;
difflickyhead=diff(lickyhead);
inicios=find(difflickyhead<0);
finales=find(difflickyhead>0);
finales=finales(2:end);
for i=2:(length(inicios)-3); %%% quitar solod debe haber un 1 en la secuencia !!!! corregir
    trials{i,1}=banderas_mod(inicios(i,1):finales(i,1),1);
    encontrar1{i,1}=find(trials{i,1}==1); 
    timemodif{i,1}=volts((inicios(i,1):finales(i,1)),4);
    
   
    if length(encontrar1{i,1})>1
     trials{i,1}=trials{i,1}(encontrar1{i,1}(1,1):(encontrar1{i,1}(2,1)),1);  
    timemodif{i,1}=timemodif{i,1}(encontrar1{i,1}(1,1):(encontrar1{i,1}(2,1)),1);
     
    end
    timealigned{i,1}=diff(timemodif{i,1});timealigned2{i,1}= vertcat(0,timealigned{i,1});
    pressizq{i,1}=trials{i,1}==4; pressizq_index{i,1}=find(pressizq{i,1}>0);timepressizq{i,1}=timealigned{i,1}((pressizq_index{i,1}-1),1); timepressizqfiber{i,1}=timemodif{i,1}((pressizq_index{i,1}-1),1);
    pressderec{i,1}=trials{i,1}==3;pressderec_index{i,1}=find(pressderec{i,1}>0);timepressdere{i,1}=timealigned{i,1}((pressderec_index{i,1}-1),1);
    recompensa{i,1}=trials{i,1}==5;recompensasum(i,:)=sum(recompensa{i,1});recompensasumindex=find(recompensasum>0);
    recompG{i,1}=trials{i,1}==6;recompGsum(i,1)=sum(recompG{i,1});
    recompNG{i,1}=trials{i,1}==7;recompNGsum(i,1)=sum(recompNG{i,1});
    castigoP1G2{i,1}=trials{i,1}==8;casigoP1Gsum(i,1)=sum(castigoP1G2{i,1});
    %castigoP2G2{i,1}=trials{i,1}==9;
    castigoP1NG2{i,1}=trials{i,1}==10;castigoP1SP(i,1)=sum(castigoP1NG2{i,1});% check that is working punish
    %casstigoP2NG2{i,1}=trials{i,1}==11;
    timeencasi{i,1}=tiempo_acum(inicios(i,1):finales(i,1),1);
    sumpresP1(i,1)=sum(pressizq{i,1});numerop1_1(i,1)=(sumpresP1(i,1)==1);numerop1_2(i,1)=(sumpresP1(i,1)==2);numerop1_3(i,1)=(sumpresP1(i,1)==3);numerop1_4(i,1)=(sumpresP1(i,1)==4);numerop1_5(i,1)=(sumpresP1(i,1)==5);numerop1_more(i,1)=(sumpresP1(i,1)>6); 
    sumpresP2(i,1)=sum(pressderec{i,1});numerop2_1(i,1)=(sumpresP2(i,1)<=2);numerop2more(i,1)=(sumpresP2(i,1)>=3);
    sum_tot(i,1)=sumpresP1(i,1)+sumpresP2(i,1);totaltrial=find(sum_tot);
    diferentialtrials{i,1}=diff(trials{i,1});
    headentry44{i,1}=trials{i,1}==1; 
    if (length(pressizq_index{i,1}))>1
    time4{i,1}=timemodif{i,1}-timemodif{i,1}(pressizq_index{i,1}(1,1),1);
    Ipi1{i,1}=timepressizq{i, 1}(2,1);
 
    end
    if length(pressizq_index{i,1})>2 
       ipi2{i,1}=timepressizq{i, 1}(3,1);     
        
    end
    if length(pressizq_index{i,1})>3    
        ipi3{i,1}=timepressizq{i, 1}(4,1);
    end
     
    
    if   length(pressizq_index{i,1})>3 && length(pressderec_index{i,1})>0 ; 
        swtichtime(i,1)=timepressdere{i, 1}(1,1);  
           
     
           
           
           
    else  
      
        swtichtime(i,1)=nan;       
    end
    
    
    if sumpresP1(i,1)>3 && sumpresP2(i,1)>0 && recompensasum(i,:)==1;
        trialscorrect{i,1}=trials{i,1};
    else
        trialscorrect{i,1}=0;
    end
    
    trialsmay2(i,1)=length(trialscorrect{i,1})>1; % seleccionando ensayos con trialss mayores a 3 p en ambas palancas
    
    %phenotipes percent related sequence 1
%     if length(timepressizq{i,1})
%         
%     end 
%     
    
    
end
recompGsum_index=find(recompGsum>0);
recompNGsum_index=find(recompNGsum>0);




trialscorrect2=find(trialsmay2>0);
trialsincorrect=find(trialsmay2==0);
totalguidedcorrect=sum(recompGsum);
totalnonguidedcorrect=sum(recompNGsum);

medianswtich=nanmedian(swtichtime);
%% 
%phenotipes percent related sequence 1

error_notcompletesequence= numerop1_1 + numerop1_2+numerop1_3;error_notcompletesequenceindex=find(error_notcompletesequence>0);
normalsequence=numerop1_4+numerop1_5;normalsequenceindex=find(normalsequence>0);
excesive_pressing=find(numerop1_more>0 & recompensasum>0 );  
sumphenotypes=length(excesive_pressing) + length (error_notcompletesequenceindex) + length(normalsequenceindex);
%percent
execivepressing1= (length(excesive_pressing))/sumphenotypes;
impulsivereact=(length (error_notcompletesequenceindex))/sumphenotypes;
normalsequence=(length(normalsequenceindex))/sumphenotypes;
%phenotipes realted lever 2 completion signal
excesivepress2=find(numerop2more>0);
normalcompletion=find(numerop2_1>0 & recompensasum>0 & sumpresP1>0);
noncompletion=find(sumpresP2>0 & recompensasum>0 & sumpresP1<1);%this is not the trial is one below
errorstattwronglever=find(sumpresP2>0 & recompensasum<1 & error_notcompletesequence<1);

error_notcompletesequenceindex;
errornoncompletemas1=error_notcompletesequenceindex+1;
errornoncompletemenos1=error_notcompletesequenceindex-1;


recompNGsumindex=find(recompNGsum>0);
recompGsumindex=find(recompGsum>0);blocks=recompNGsumindex(1,1):5:length(recompGsumindex);blocks=blocks';

for v=1:length(error_notcompletesequenceindex);
  errorNG{v,1}= find(recompNGsumindex==errornoncompletemas1(v,1));
  errorNG2(v,1)=length(errorNG{v,1});
  
end
errorNG3nc=find(errorNG2>0);
errorG=find(errorNG2<1);
PerformanceG=length(recompGsumindex)/(length(recompGsumindex)+length(errorG));
PerformanceNG=length(recompNGsumindex)/((length(recompNGsumindex)+length(errorNG3nc))+length(errorstattwronglever));

sumphenotipespress2=(length(excesivepress2))+ (length(noncompletion)) + (length(normalcompletion))+(length(errorstattwronglever));
Excesivecompletion=(length(excesivepress2))/sumphenotipespress2;
normalcompletion2=length(normalcompletion)/sumphenotipespress2;;
impulsivecompletion1=length(noncompletion)/sumphenotipespress2;;

numberrewardsindex=find(recompensasum>0);
totaltrials=(length(numberrewardsindex))+ (length(error_notcompletesequenceindex)) +(length(errorstattwronglever));

performance=(length(numberrewardsindex))/totaltrials;
errorincompletesequence=(length(error_notcompletesequenceindex))/totaltrials;
errorstartwronlever2=(length(errorstattwronglever))/totaltrials;





%%
figure(1)
hold on
plot(volts(:,4),press1,'b+');
plot(volts(:,4),press2,'r+');
plot(volts(:,4),reward,'b.');
plot(volts(:,4),headentry,'*g');
trialscorrected=find(trialsmay2>0);



%% performance
datas3{1,1}=performance;
datas3{1,2}=errorincompletesequence;
datas3{1,3}=errorstartwronlever2;

%% phenotipes lever1
datas3{1,4}=execivepressing1;
datas3{1,5}=impulsivereact;
datas3{1,6}=normalsequence;

%% phenotipes lever2
datas3{1,7}=Excesivecompletion;
datas3{1,8}=impulsivecompletion1;
datas3{1,9}=normalcompletion2;

datas3{1,10}=PerformanceG;
datas3{1,11}=PerformanceNG;
datas3{1,12}=medianswtich;
datas3{2,2}=timespress2;
datas3{2,3}=timespressdiff;
datas3{3,3}=timespress2+timealinfotometry;
%datas3{2,5}=
datas3{2,1}=timeleverP1bb
datas3{3,1}=timespress1diff;
datas3{4,1}=timespress1diff+timealinfotometry;





%% Figure 2 completed trials

if length(normalsequenceindex>2)

figure (2)
hold on
for pp=1:length(normalsequenceindex);
     if length(time4{normalsequenceindex(pp,1),1})>0
    Latencycorrect(pp,1)= abs(time4{normalsequenceindex(pp,1),1}(1,1));
    Switchcorrect(pp,1)=swtichtime(normalsequenceindex(pp,1),1);
    plot((time4{normalsequenceindex(pp,1),1}),(pp*pressizq{normalsequenceindex(pp,1),1}), 'b+');
    plot((time4{normalsequenceindex(pp,1),1}),(pp*pressderec{normalsequenceindex(pp,1),1}), 'r+');
    plot((time4{normalsequenceindex(pp,1),1}),(pp*recompensa{normalsequenceindex(pp,1),1}), 'gv')
    plot((time4{normalsequenceindex(pp,1),1}),(pp*recompG{normalsequenceindex(pp,1),1}), 'cv')
    plot((time4{normalsequenceindex(pp,1),1}),(pp*recompNG{normalsequenceindex(pp,1),1}), 'mv')
    plot((time4{normalsequenceindex(pp,1),1}),(pp*headentry44{normalsequenceindex(pp,1),1}), 'b.')
    
    
    %plot((time4{trialscorrect2(pp,1),1}),(pp*time4{trialscorrect2(pp,1)}(1,1),1}), 'm.')
    axis( [-10 10 1 120])
    xlabel('time (sec)');
    ylabel('sequences');
    title('trials completed')
    ipi1correct(pp,1)=(Ipi1{normalsequenceindex(pp,1),1});
    ipi2correct(pp,1)=(ipi2{normalsequenceindex(pp,1),1});
    ipi3correct(pp,1)=(ipi3{normalsequenceindex(pp,1),1});
    startcorrect(pp,1)=timemodif{normalsequenceindex(pp,1),1}(pressizq_index{normalsequenceindex(pp,1),1}(1,1),1);
    secondcorrect(pp,1)=timemodif{normalsequenceindex(pp,1),1}(pressizq_index{normalsequenceindex(pp,1),1}(2,1),1);
    endcorrect(pp,1)=timemodif{normalsequenceindex(pp,1),1}(pressizq_index{normalsequenceindex(pp,1),1}(end,1),1);
    if sum (pressderec{normalsequenceindex(pp,1),1})>0
         completionlever(pp,1)= timemodif{normalsequenceindex(pp,1),1}(pressderec_index {normalsequenceindex(pp,1),1}(1,1),1);
    else
         completionlever(pp,1)=nan;
        
    end
    
    
    if  recompGsum(normalsequenceindex(pp,1),1)> 0
       
    ipi1correctG(pp,1)=(Ipi1{normalsequenceindex(pp,1),1});
    ipi2correctG(pp,1)=(ipi2{normalsequenceindex(pp,1),1});
    ipi3correctG(pp,1)=(ipi3{normalsequenceindex(pp,1),1});
    startcorrectG(pp,1)=timemodif{normalsequenceindex(pp,1),1}(pressizq_index{normalsequenceindex(pp,1),1}(1,1),1);
    secondcorrectG(pp,1)=timemodif{normalsequenceindex(pp,1),1}(pressizq_index{normalsequenceindex(pp,1),1}(2,1),1);
    endcorrectG(pp,1)=timemodif{normalsequenceindex(pp,1),1}(pressizq_index{normalsequenceindex(pp,1),1}(end,1),1);
    completionleverG(pp,1)= completionlever(pp,1);
   
    completionleverNG(pp,1)=nan;
        
    else 
    ipi1correctSP(pp,1)=(Ipi1{normalsequenceindex(pp,1),1});
    ipi2correctSP(pp,1)=(ipi2{normalsequenceindex(pp,1),1});
    ipi3correctSP(pp,1)=(ipi3{normalsequenceindex(pp,1),1});
    startcorrectSP(pp,1)=timemodif{normalsequenceindex(pp,1),1}(pressizq_index{normalsequenceindex(pp,1),1}(1,1),1);
    secondcorrectSP(pp,1)=timemodif{normalsequenceindex(pp,1),1}(pressizq_index{normalsequenceindex(pp,1),1}(2,1),1);
    endcorrectSP(pp,1)=timemodif{normalsequenceindex(pp,1),1}(pressizq_index{normalsequenceindex(pp,1),1}(end,1),1);
    completionleverNG(pp,1)= completionlever(pp,1);
    completionleverG(pp,1)=nan;
        
    end
    
    
     end
end 

datas3{2,5}=startcorrect;
datas3{3,5}=startcorrect+timealinfotometry;
datas3{4,5}=secondcorrect;
datas3{5,5}=secondcorrect+timealinfotometry;
datas3{6,5}=endcorrect;
datas3{7,5}=endcorrect+timealinfotometry;
datas3{8,5}=completionlever;
datas3{9,5}=completionlever+timealinfotometry;

%Guided chain
datas3{10,5}=startcorrectG
datas3{11,5}=startcorrectG+timealinfotometry;
datas3{12,5}=secondcorrectG
datas3{13,5}=secondcorrectG+timealinfotometry;
datas3{14,5}=endcorrectG
datas3{15,5}=endcorrectG+timealinfotometry;
datas3{16,5}=completionleverG
datas3{17,5}=completionleverG+timealinfotometry;
%selfpace chain
datas3{18,5}=startcorrectSP
datas3{19,5}=startcorrectSP+timealinfotometry;
datas3{20,5}=secondcorrectSP
datas3{21,5}=secondcorrectSP+timealinfotometry;
datas3{22,5}=endcorrectSP
datas3{23,5}=endcorrectSP+timealinfotometry;
datas3{24,5}=completionleverNG
datas3{25,5}=completionleverNG+timealinfotometry;


datas3{1,13}=nanmedian(ipi1correct);
datas3{1,16}=nanmedian(ipi2correct);
datas3{1,19}=nanmedian(ipi3correct);
datas3{1,21}=nanmedian(Latencycorrect);
datas3{1,24}=nanmedian(Switchcorrect);

 
[Where] = Where_file(filename);
SaveName= [strcat(Where,'figure_completed_trials',name)];
     savefig(SaveName);
end
%% figure 3 early swtich  impulsivity
if length(error_notcompletesequenceindex>2)
figure (3)
hold on
for pp=1:length(error_notcompletesequenceindex);
    if error_notcompletesequenceindex(pp,1)<length(time4)
    if length(time4{error_notcompletesequenceindex(pp,1),1})>0
     Latencyimpulsiv(pp,1)= abs(time4{error_notcompletesequenceindex(pp,1),1}(1,1));
     Switchimpul(pp,1)= time4{error_notcompletesequenceindex(pp,1),1}(end,1)-time4{ error_notcompletesequenceindex(pp,1),1}(pressizq_index{error_notcompletesequenceindex(pp,1),1}(end,1));
     
    plot((time4{error_notcompletesequenceindex(pp,1),1}),(pp*pressizq{error_notcompletesequenceindex(pp,1),1}), 'b+');
    plot((time4{error_notcompletesequenceindex(pp,1),1}),(pp*pressderec{error_notcompletesequenceindex(pp,1),1}), 'r+');
    plot((time4{error_notcompletesequenceindex(pp,1),1}),(pp*recompensa{error_notcompletesequenceindex(pp,1),1}), 'gv')
    plot((time4{error_notcompletesequenceindex(pp,1),1}),(pp*recompG{error_notcompletesequenceindex(pp,1),1}), 'cv')
    plot((time4{error_notcompletesequenceindex(pp,1),1}),(pp*recompNG{error_notcompletesequenceindex(pp,1),1}), 'mv')
    plot((time4{error_notcompletesequenceindex(pp,1),1}),(pp*headentry44{error_notcompletesequenceindex(pp,1),1}), 'b.')
  
    %plot((time4{trialscorrect2(pp,1),1}),(pp*time4{trialscorrect2(pp,1)}(1,1),1}), 'm.')
    axis( [-10 10 1 120])
    xlabel('time (sec)');
    ylabel('sequences');
    title('impulsive trials')
   
        
     if (length(Ipi1{error_notcompletesequenceindex(pp,1),1})>0) 
     ipi1impul(pp,1)=Ipi1{error_notcompletesequenceindex(pp,1),1};
      else  ipi1impul(pp,1)=nan; end
   if (length(ipi2{error_notcompletesequenceindex(pp,1),1})>0); ipi2iimpul(pp,1)=ipi2{error_notcompletesequenceindex(pp,1),1};else ipi2iimpul(pp,1)=nan;end 
    %if (length(ipi3{error_notcompletesequenceindex(pp,1),1})>0);ipi3impul(pp,1)=(ipi3{error_notcompletesequenceindex(pp,1),1});else ipi3impul(pp,1)=nan;end
          startimpuls(pp,1)=timemodif{error_notcompletesequenceindex(pp,1),1}(pressizq_index{error_notcompletesequenceindex(pp,1),1}(1,1),1);
         secondimpuls(pp,1)=timemodif{error_notcompletesequenceindex(pp,1),1}(pressizq_index{error_notcompletesequenceindex(pp,1),1}(2,1),1);
          endimpuls(pp,1)=timemodif{error_notcompletesequenceindex(pp,1),1}(pressizq_index{error_notcompletesequenceindex(pp,1),1}(end,1),1);

    end
    end
end 

datas3{2,6}=startimpuls;
datas3{3,6}=startimpuls+timealinfotometry;
datas3{4,6}=secondimpuls;
datas3{5,6}=secondimpuls+timealinfotometry;
datas3{6,6}=endimpuls;
datas3{7,6}=endimpuls+timealinfotometry;


datas3{1,14}=nanmedian(ipi1impul);
datas3{1,22}=nanmedian(Latencyimpulsiv);
datas3{1,17}=nanmedian(ipi2iimpul);
datas3{1,25}=nanmedian(Switchimpul);
[Where] = Where_file(filename);
SaveName= [strcat(Where,'figure_impulsive_trials',name)];
     savefig(SaveName);
end
%% overpressing    compulsivity

if length(excesive_pressing)>2

figure (4)
hold on
for pp=1:length(excesive_pressing);
     if length(time4{excesive_pressing(pp,1),1})>0
      Latencycompul(pp,1)= abs(time4{excesive_pressing(pp,1),1}(1,1));
     Switchcompul(pp,1)=swtichtime(excesive_pressing(pp,1),1);
        
         
    plot((time4{excesive_pressing(pp,1),1}),(pp*pressizq{excesive_pressing(pp,1),1}), 'b+');
    plot((time4{excesive_pressing(pp,1),1}),(pp*pressderec{excesive_pressing(pp,1),1}), 'r+');
    plot((time4{excesive_pressing(pp,1),1}),(pp*recompensa{excesive_pressing(pp,1),1}), 'gv')
    plot((time4{excesive_pressing(pp,1),1}),(pp*recompG{excesive_pressing(pp,1),1}), 'cv')
    plot((time4{excesive_pressing(pp,1),1}),(pp*recompNG{excesive_pressing(pp,1),1}), 'mv')
    plot((time4{excesive_pressing(pp,1),1}),(pp*headentry44{excesive_pressing(pp,1),1}), 'b.')
    %plot((time4{trialscorrect2(pp,1),1}),(pp*time4{trialscorrect2(pp,1)}(1,1),1}), 'm.')
    axis( [-10 10 1 120])
    xlabel('time (sec)');
    ylabel('sequences');
    title('compulsive trials');
    if (length(Ipi1{excesive_pressing(pp,1),1})>0);
    ipi1compul(pp,1)=Ipi1{excesive_pressing(pp,1),1};
    else  ipi1compul(pp,1)=nan; end
    if (length(ipi2{excesive_pressing(pp,1),1})>0); ipi2icomp(pp,1)=ipi2{excesive_pressing(pp,1),1};else ipi1compul(pp,1)=nan;end 
    if (length(ipi3{excesive_pressing(pp,1),1})>0);ipi3imcompuls(pp,1)=(ipi3{excesive_pressing(pp,1),1});else ipi1compul(pp,1)=nan;end
     
  startoverpress(pp,1)=timemodif{excesive_pressing(pp,1),1}(pressizq_index{excesive_pressing(pp,1),1}(1,1),1);
  secondoverpress(pp,1)=timemodif{excesive_pressing(pp,1),1}(pressizq_index{excesive_pressing(pp,1),1}(2,1),1);
  endoverpress(pp,1)=timemodif{excesive_pressing(pp,1),1}(pressizq_index{excesive_pressing(pp,1),1}(end,1),1);
   if sum (pressderec{excesive_pressing(pp,1),1})>0
         completionleverComp(pp,1)= timemodif{excesive_pressing(pp,1),1}(pressderec_index {excesive_pressing(pp,1),1}(1,1),1);
    else
         completionleverComp(pp,1)=nan;
        
    end
  
     end
end 



[Where] = Where_file(filename);
SaveName= [strcat(Where,'figure_compulsive_trials',name)];
     savefig(SaveName);

     datas3{2,7}=startoverpress;
     datas3{3,7}=startoverpress+timealinfotometry;
     datas3{4,7}=secondoverpress;
     datas3{5,7}=secondoverpress+timealinfotometry;
     datas3{6,7}=endoverpress;
     datas3{7,7}=endoverpress+timealinfotometry;
     datas3{8,7}=completionleverComp;
     datas3{9,7}=completionleverComp + timealinfotometry;
     
     
     datas3{1,15}=nanmedian(ipi1compul);
     datas3{1,18}=nanmedian(ipi2icomp);
     datas3{1,20}=nanmedian(ipi3imcompuls);
     datas3{1,23}=nanmedian(Latencycompul);
     datas3{1,26}=nanmedian(Switchcompul);
end


    
   
  
     
     
end 