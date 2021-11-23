close all;
clear all;

%% Extract photometry data from TDT Data Tank

%path_to_dropbox = '/Users/elizabethholly/Desktop/PhotometryData/SSC310_1023'; % point to folder where data is
%tankdir = path_to_dropbox;
%tankname = 'SSC310-171018-113156';
%blockname = 'SSC310_1023';
%new_direc = strcat(tankdir, '\',tankname, '\', blockname);
%cd (new_direc);
%storenames = {'470B' '405B' 'Ms1/' 'Ms1\' 'Ms2/' 'Ms2\' };



    filepath='C:\Users\edgar\Desktop\Nueva carpeta (3)\5_1_and_6_1-210608-155017\5_1_and_6_1';
    %exptname=base; filename=dir+exptname; savename=filename+100Hz
%     exptnameDef=strsplit(filepath,'/');
%     exptnameDef=char(exptnameDef(length(exptnameDef)));  %the savename becomes the folder name
%     exptnameDef = {exptnameDef}; %default expt name saved as a cell array of strings
%     exptname = inputdlg('experiment name?','Input',1,exptnameDef,'on'); exptname = exptname{1};
%     mkdir(exptname); addpath(exptname);% make a folder name for the m file output
%     filename=strcat(exptname,'/',exptname); %savename dir becomes destination and filename



filename = 'SSC310_1023.mat'; %filename where entire workspace can be saved
figname = 'SSC310_1023signal.fig';  %figname for three signal graphs
figname1 = 'SSC310_1023final.fig';  %figname for final data graph
storenames = {'40AA' '45AA'  'Ms4\' 'Ms4/' }; 
 
% extract
for k = 1:numel(storenames)
  storename = storenames{k};
  %S{k} = tdt2mat2(tankdir, tankname, blockname, storename);
   S{k} = TDT_MATconvertBDH(filepath, storename); %BHedit
end



%% Massage data and get time stamps
% Get LMag data as a vector (repeat for each channel)
dat1 = S{1};
dat1.data = reshape(dat1.data', [],1); % unwrap data from m x 256 array
dat1.data = dat1.data(1:end-5000);
dat2 = S{2};
dat2.data = reshape(dat2.data', [],1); % unwrap data from m x 256 array
dat2.data = dat2.data(1:end-5000);

% Get LMag timestamps (use chani1 - timestamps should be the same for all LMag channels
ts = dat1.timestamps;
t_rec_start = ts(1);
samplingRate = dat1.sampling_rate; % TDT returns 'true' sample rate.

ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:dat1.npoints-1)*(1./dat1.sampling_rate));
ts = reshape(ts',[],1);
ts = ts(1:end-5000);

% Get TTL input timestamps
TTL1 = S{3}.timestamps - t_rec_start;
TTL2 = S{4}.timestamps - t_rec_start;
% TTL3 = S{5}.timestamps - t_rec_start;
% TTL4 = S{6}.timestamps - t_rec_start;

%% Single fit - if both channels change consistently for the whole session

% Smooth dat1 and dat2, fit dat2 (control) to dat1 (signal)

dat1 = filtfilt(ones(1,100)/100,1, dat1.data);
dat2 = filtfilt(ones(1,100)/100,1, dat2.data);

tsOrig=ts;
[dat1,ts]=resample(dat1,tsOrig,40,10,10);
dat2=resample(dat2,tsOrig,40,10,10);
samplingRate = 40;

debleach_flag=1;
[dfof, mod_dat1, fit_dat1, offset_dat1]=debleachBDH(ts',dat1,debleach_flag);
[dfof_control, mod_dat2, fit_dat2, offset_dat2] =debleachBDH(ts',dat2,debleach_flag);

 f(1) = figure('Color',[1 1 1]);
 zoom on;
 subplot (2, 1, 1);
 plot(ts/60,dat1);
 hold on
 plot(ts/60,fit_dat1+offset_dat1,'r');
 plot(ts/60,mod_dat1,'g');
 legend('Signal','Exponential Fit','Debleached');
 legend BOXOFF;
 xlabel('Time(s)');
 ylabel('Fluorescence');
 title('Debleaching Results for 470');
 
%  f(2) = figure('Color',[1 1 1]);
%  zoom on;
%  subplot (2, 1, 2);
 %plot(ts/60,dat2);
 %hold on
 %plot(ts/60,fit_dat2+offset_dat2,'r');
 %plot(ts/60,mod_dat2,'g');
 %legend('Signal','Exponential Fit','Debleached');
 %legend BOXOFF;
 %xlabel('Time(s)');
 %ylabel('Fluorescence');
 %title('Debleaching Results for 405');

dfofCorr=subtract_refBDH(ts',dfof,dfof_control,'None'); %for some stupid reason, it still says subtract, but it isn't. on the other one, it is.
dfofCorr_sub=subtract_refBDH(ts',dfof,dfof_control,'Subtract');

[controlFit] = controlFit (dat1, dat2);
%from here on, what ope did was take three different ways to calculate
%dfof, there was the one from ES originally, that is the normdat from her
%code. anything with just normdat with no db is that. anything with db is
%debleached. z means it got z scored, with no moving window. modz means
%moving window. anythign with diff is a differenced time series. 
% Get delta F/F using controlFit
[normDat] = deltaFF (dat1, controlFit);

%[normDat] = deltaFF (dat1, controlFit);
normDat_DB = 100*dfofCorr';
normDat_DBZ = (normDat_DB-mean(normDat_DB))/std(normDat_DB);
% normDat_DBZ_diff = diff(normDat_DBZ);

%[normDat] = deltaFF (dat1, controlFit);
normDat_DB_sub = 100*dfofCorr_sub';
normDat_DBZ_sub = (normDat_DB_sub-mean(normDat_DB_sub))/std(normDat_DB_sub);
% normDat_DBZ_sub_diff = diff(normDat_DBZ_sub);
% 
%calculate local means and std in 10s windows
% moving_mean = zeros(1,length(normDat));
% moving_sd = zeros(1,length(normDat));
% moving_mean(1:10*samplingRate) = mean(normDat(1:10*samplingRate));
% moving_sd(1:10*samplingRate) = std(normDat(1:10*samplingRate));
% moving_mean(end-(10*samplingRate):end) = mean(normDat(end-(10*samplingRate):end));
% moving_sd(end-(10*samplingRate):end) = std(normDat(end-(10*samplingRate):end));
% for i = (10*samplingRate)+1:(length(normDat)-(10*samplingRate))
%     moving_mean(i) = mean(normDat((i-10*samplingRate):(i+10*samplingRate)));
%     moving_sd(i) = std(normDat((i-10*samplingRate):(i+10*samplingRate)));
% end
% 
% %z-score data using local means and standard deviations
% normDat_modZ = zeros(1,length(normDat_DB));
% for i = 1:length(normDat)
%     normDat_modZ(i) = (normDat(i)-moving_mean(i))/moving_sd(i);
% end
% 
% normDat_modZ = normDat_modZ';

%calculate local means and std in 10s windows for debleached data
%with no subtraction
% moving_mean = zeros(1,length(normDat_DB));
% moving_sd = zeros(1,length(normDat_DB));
% moving_mean(1:10*samplingRate) = mean(normDat_DB(1:10*samplingRate));
% moving_sd(1:10*samplingRate) = std(normDat_DB(1:10*samplingRate));
% moving_mean(end-(10*samplingRate):end) = mean(normDat_DB(end-(10*samplingRate):end));
% moving_sd(end-(10*samplingRate):end) = std(normDat_DB(end-(10*samplingRate):end));
% for i = (10*samplingRate)+1:(length(normDat_DB)-(10*samplingRate))
%     moving_mean(i) = mean(normDat_DB((i-10*samplingRate):(i+10*samplingRate)));
%     moving_sd(i) = std(normDat_DB((i-10*samplingRate):(i+10*samplingRate)));
% end
% 
% %z-score data using local means and standard deviations
% normDat_DB_modZ = zeros(1,length(normDat_DB));
% for i = 1:length(normDat_DB)
%     normDat_DB_modZ(i) = (normDat_DB(i)-moving_mean(i))/moving_sd(i);
% end
% 
% normDat_DB_modZ = normDat_DB_modZ';
% % 
% %calculate local means and std in 10s  windows for debleached data with
% %subtraction of 405
% moving_mean = zeros(1,length(normDat_DB_sub));
% moving_sd = zeros(1,length(normDat_DB_sub));
% moving_mean(1:10*samplingRate) = mean(normDat_DB_sub(1:10*samplingRate));
% moving_sd(1:10*samplingRate) = std(normDat_DB_sub(1:10*samplingRate));
% moving_mean(end-(10*samplingRate):end) = mean(normDat_DB_sub(end-(10*samplingRate):end));
% moving_sd(end-(10*samplingRate):end) = std(normDat_DB_sub(end-(10*samplingRate):end));
% for i = (10*samplingRate)+1:(length(normDat_DB_sub)-(10*samplingRate))
%     moving_mean(i) = mean(normDat_DB_sub((i-10*samplingRate):(i+10*samplingRate)));
%     moving_sd(i) = std(normDat_DB_sub((i-10*samplingRate):(i+10*samplingRate)));
% end
% 
% %z-score data using local means and standard deviations
% normDat_DB_sub_modZ = zeros(1,length(normDat_DB_sub));
% for i = 1:length(normDat_DB_sub)
%     normDat_DB_sub_modZ(i) = (normDat_DB_sub(i)-moving_mean(i))/moving_sd(i);
% end
% 
% normDat_DB_sub_modZ = normDat_DB_sub_modZ';

%this is old. it will plot everything, thats just so that you can inspect
%it to see what looks good. what you can tell a lot is whether the
%subtracting 405 adds noise, whether you want to, etc. 
% f(3) = figure('Color',[1 1 1]);
% zoom on;
% subplot (2, 1, 1);
% plot(ts(51:end)/60, normDat(51:end));
% xlabel('Time(s)');
% ylabel('\DeltaF/F (%)');
% title('\DeltaF/F - No Transformations');

% f(4) = figure('Color',[1 1 1]);
% zoom on;
% subplot (2, 1, 2);
% plot(ts(51:end)/60, normDat_modZ(51:end));
% xlabel('Time(s)');
% ylabel('\DeltaF/F (z-score)');
% title('\DeltaF/F Z-Score (Local Mean/Standard Deviation); Original Data)');

f(5) = figure('Color',[1 1 1]);
zoom on;
subplot (4, 1, 1);
plot(ts/60, normDat_DB);
xlabel('Time(min)');
ylabel('\DeltaF/F (%)');
title('\DeltaF/F Debleached');

f(6) = figure('Color',[1 1 1]);
zoom on;
subplot (4, 1, 2);
plot(ts/60, normDat_DBZ);
xlabel('Time(min)');
ylabel('\DeltaF/F (z-score)');
title('\DeltaF/F Z-Score (Global Standard Deviation); Debleached Data)');

%  f(7) = figure('Color',[1 1 1]);
%  zoom on;
%  subplot (4, 1, 3);
%  plot(ts/60, normDat_DB_modZ);
%  xlabel('Time(s)');
%  ylabel('\DeltaF/F (z-score)');
%  title('\DeltaF/F Z-Score (Local Mean/Standard Deviation); Debleached Data)');

% f(8) = figure('Color',[1 1 1]);
% zoom on;
% subplot (4, 1, 4);
% plot(ts(2:end)/60, normDat_DBZ_diff);
% xlabel('Time(s)');
% ylabel('\DeltaF/F (z-score)');
% title('Differencing of Debleached Z-Score Data)');

% f(9) = figure('Color',[1 1 1]);
% zoom on;
% %subplot (4, 1, 1);
% plot(ts/60, normDat_DB_sub);
% xlabel('Time(s)');
% ylabel('\DeltaF/F (%)');
% title('\DeltaF/F Debleached w/ Subtraction of 405');

% f(10) = figure('Color',[1 1 1]);
% zoom on;
% subplot (4, 1, 2);
% plot(ts/60, normDat_DBZ_sub);
% xlabel('Time(s)');
% ylabel('\DeltaF/F (z-score)');
% title('\DeltaF/F Z-Score (Global Standard Deviation); Debleached Data w/ 405 Subtraction');

% f(11) = figure('Color',[1 1 1]);
% zoom on;
% subplot (4, 1, 3);
% plot(ts/60, normDat_DB_sub_modZ);
% xlabel('Time(s)');
% ylabel('\DeltaF/F (z-score)');
% title('\DeltaF/F Z-Score (Local Mean/Standard Deviation); Debleached Data w/ 405 Subtraction');

% f(12) = figure('Color',[1 1 1]);
% zoom on;
% subplot (4, 1, 4);
% plot(ts(2:end)/60, normDat_DBZ_sub_diff);
% xlabel('Time(s)');
% ylabel('\DeltaF/F (z-score)');
% title('Differencing of Debleached Z-Score Data w/ 405 Subtraction');


save('SSC312_1024_peak_analysis.mat');
%% Save raw, fit and dFF figs to single file

% savefig (gcf, figname);

% %this is the peakfinder stuff
% [DBZpeakloc,DBZrealpeakMag, DBZrealpeakwidth, DBZrealpeakprom, DBZrealpeakheight, f(5)] = peakfinder_EHedit_10s_withmodeBL(normDat_DBZ, samplingRate, ts);
% zoom on;
% DBZrealpeakMag=transpose(DBZrealpeakMag); DBZrealpeakprom=transpose(DBZrealpeakprom); DBZrealpeakwidth=transpose(DBZrealpeakwidth); DBZrealpeakheight=transpose(DBZrealpeakheight);
% DBZpeakinfo=[DBZpeakloc, DBZrealpeakMag, DBZrealpeakprom, DBZrealpeakwidth, DBZrealpeakheight];
% 
% [DBZ_subpeakloc,DBZ_subrealpeakMag, DBZ_subrealpeakwidth, DBZ_subrealpeakprom, DBZ_subrealpeakheight, f(5)] = peakfinder_EHedit_10s_withmodeBL(normDat_DBZ_sub, samplingRate, ts);
% zoom on;
% DBZ_subrealpeakMag=transpose(DBZ_subrealpeakMag); DBZ_subrealpeakprom=transpose(DBZ_subrealpeakprom); DBZ_subrealpeakwidth=transpose(DBZ_subrealpeakwidth); DBZ_subrealpeakheight=transpose(DBZ_subrealpeakheight);
% DBZ_subpeakinfo=[DBZ_subpeakloc, DBZ_subrealpeakMag, DBZ_subrealpeakprom, DBZ_subrealpeakwidth, DBZ_subrealpeakheight];
% 
% [DB_subpeakloc,DB_subrealpeakMag, DB_subrealpeakwidth, DB_subrealpeakprom, DB_subrealpeakheight, f(5)] = peakfinder_EHedit_10s_withmodeBL(normDat_DB_sub, samplingRate, ts);
% zoom on;
% DB_subrealpeakMag=transpose(DB_subrealpeakMag); DB_subrealpeakprom=transpose(DB_subrealpeakprom); DB_subrealpeakwidth=transpose(DB_subrealpeakwidth); DB_subrealpeakheight=transpose(DB_subrealpeakheight);
% DB_subpeakinfo=[DB_subpeakloc, DB_subrealpeakMag, DB_subrealpeakprom, DB_subrealpeakwidth, DB_subrealpeakheight];
% 
% [DBpeakloc,DBrealpeakMag, DBrealpeakwidth, DBrealpeakprom, DBrealpeakheight, f(5)] = peakfinder_EHedit_10s_withmodeBL(normDat_DB, samplingRate, ts);
% zoom on;
% DBrealpeakMag=transpose(DBrealpeakMag); DBrealpeakprom=transpose(DBrealpeakprom); DBrealpeakwidth=transpose(DBrealpeakwidth); DBrealpeakheight=transpose(DBrealpeakheight);
% DBpeakinfo=[DBpeakloc, DBrealpeakMag, DBrealpeakprom, DBrealpeakwidth, DBrealpeakheight];

%% make  arrays based on behavior timestamps, then average rows into a mean vector for plotting
nSecPrev = 5; %change to make different window 
nSecPost = 20; 
% convert seconds to TDT timestamps
nTsPrev = round (nSecPrev * samplingRate);
nTsPost = round (nSecPost * samplingRate);
%%%%%%%%%
%% make  arrays for behavior timestamps 


% make PSTH for TTL1
nTTL1 = size(TTL1,1);
PsthArray_TTL1 = NaN(nTTL1,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nTTL1
    thisTime = TTL1(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_TTL1(i,:) = processPhotDataRow_normDat(normDat_DB, thisIndex, nTsPrev, nTsPost);
end
err_TTL1 = (nanstd(PsthArray_TTL1))/sqrt(size(PsthArray_TTL1,1));
Psth_TTL1 = nanmean(PsthArray_TTL1);


% make PSTH for TTL2
nTTL2 = size(TTL2,1);
PsthArray_TTL2 = NaN(nTTL2,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nTTL2
    thisTime = TTL2(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_TTL2(i,:) = processPhotDataRow_normDat(normDat_DB, thisIndex, nTsPrev, nTsPost);
end
err_TTL2 = (nanstd(PsthArray_TTL2))/sqrt(size(PsthArray_TTL2,1));
Psth_TTL2 = nanmean(PsthArray_TTL2);



% % make PSTH for TTL3
% nTTL3 = size(TTL3,1);
% PsthArray_TTL3 = NaN(nTTL3,nTsPrev+nTsPost+1); % preallocate arrays for speed
% for i = 1:nTTL3
%     thisTime = TTL3(i);
%     thisIndex = round((thisTime*samplingRate))+1;
%     PsthArray_TTL3(i,:) = processPhotDataRow_normDat(normDat_DB_sub, thisIndex, nTsPrev, nTsPost);
% end
% err_TTL3 = (nanstd(PsthArray_TTL3))/sqrt(size(PsthArray_TTL3,1));
% Psth_TTL3 = nanmean(PsthArray_TTL3);
% 
% % make PSTH for TTL4
% nTTL4 = size(TTL4,1);
% PsthArray_TTL4 = NaN(nTTL4,nTsPrev+nTsPost+1); % preallocate arrays for speed
% for i = 1:nTTL4
%     thisTime = TTL4(i);
%     thisIndex = round((thisTime*samplingRate))+1;
%     PsthArray_TTL4(i,:) = processPhotDataRow_normDat(normDat_DB_sub, thisIndex, nTsPrev, nTsPost);
% end
% err_TTL4 = (nanstd(PsthArray_TTL4))/sqrt(size(PsthArray_TTL4,1));
% Psth_TTL4 = nanmean(PsthArray_TTL4);

    
%% set values used for all plots if desired
ymin = -2;
ymax = 2; 

%% plot all TTLs on same graph

totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
Err1Pos = Psth_TTL1 + err_TTL1;
Err1Neg = Psth_TTL1 - err_TTL1;
Err2Pos = Psth_TTL2 + err_TTL2;
Err2Neg = Psth_TTL2 - err_TTL2;
% Err3Pos = Psth_TTL3 + err_TTL3;
% Err3Neg = Psth_TTL3 - err_TTL3;
% Err4Pos = Psth_TTL4 + err_TTL4;
% Err4Neg = Psth_TTL4 - err_TTL4;
f(13) = figure('Color',[1 1 1]);
hold on

    %o = 0.5;
    %r1 = 255;
    %g1 = 255;
    %b1 = 215;
    %rgb1_o = opacity (o, r1, g1, b1);
    %fill([timeAxis, fliplr(timeAxis)],[Err1Pos, fliplr(Err1Neg)], rgb1_o, 'EdgeColor', 'none');
    %h1= plot (timeAxis,Psth_ITI,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

h1f=fill([timeAxis, fliplr(timeAxis)],[Err1Pos, fliplr(Err1Neg)], 'g', 'EdgeColor', 'none');
set(h1f,'facealpha',.5);
h1= plot (timeAxis,Psth_TTL2,'Color','g', 'Linewidth', 3);



    %r2 = 232;
    %g2 = 155;
    %b2 = 215;
    %rgb2_o = opacity (o, r2, g2, b2);
    %fill([timeAxis, fliplr(timeAxis)],[Err2Pos, fliplr(Err2Neg)], rgb2_o, 'EdgeColor', 'none');
    %h2= plot (timeAxis,Psth_INITIATION,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);
h2f=fill([timeAxis, fliplr(timeAxis)],[Err2Pos, fliplr(Err2Neg)], 'b', 'EdgeColor', 'none');
set(h2f,'facealpha',.5);
h2= plot (timeAxis,Psth_TTL1,'Color','b', 'Linewidth', 3);


    %r3 = 155;
    %g3 = 215;
    %b3 = 232;
    %rgb3_o = opacity (o, r3, g3, b3);
    %fill([timeAxis, fliplr(timeAxis)],[Err3Pos, fliplr(Err3Neg)], rgb3_o, 'EdgeColor', 'none');
    %h3= plot (timeAxis,Psth_CORRECTPRESS,'Color',[(r3/255) (g3/255) (b3/255)], 'Linewidth', 3);
% h3f=fill([timeAxis, fliplr(timeAxis)],[Err3Pos, fliplr(Err3Neg)], 'r', 'EdgeColor', 'none');
% set(h3f,'facealpha',.5);
% h3= plot (timeAxis,Psth_TTL3,'Color','r', 'Linewidth', 3);
%     
    
    %r4 = 183;
    %g4 = 232;
    %b4 = 155;
    %rgb4_o = opacity (o, r4, g4, b4);
    %fill([timeAxis, fliplr(timeAxis)],[Err4Pos, fliplr(Err4Neg)], rgb4_o, 'EdgeColor', 'none');
    %h4= plot (timeAxis,Psth_REWARD,'Color',[(r4/255) (g4/255) (b4/255)], 'Linewidth', 3);
% h4f=fill([timeAxis, fliplr(timeAxis)],[Err4Pos, fliplr(Err4Neg)], 'm', 'EdgeColor', 'none');
% set(h4f,'facealpha',0.5)
% h4= plot (timeAxis,Psth_TTL4,'Color','m', 'Linewidth', 3);



h5= line('XData', [(-1 * nSecPrev),nSecPost], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
h6= line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')

%labels, legend, make pretty, size
xlabel('Time (s)', 'FontSize', 20);
ylabel('\DeltaF/F (%)', 'FontSize', 20);
legend ([h1 h2 ], 'TTL1', 'TTL2',  'orientation', 'vertical', 'Location', 'SouthWest');
legend BOXOFF;
xlim ([(-1 * nSecPrev),nSecPost]);
ylim ([ymin, ymax]);
set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'Fontsize', 18, ...
  'LineWidth'   , 2         );


set(gca, 'layer', 'top');
set(gcf, 'Position', [0 0 700 500])
figure(20)
image(PsthArray_TTL1,'CDataMapping','scaled');
figure(21)
image(PsthArray_TTL2,'CDataMapping','scaled');
%% Save workspace variables

savefig (f, figname1, 'compact');
clear h;
clear f;
save (filename);






