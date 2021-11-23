function [realpeakloc,realpeakMag, realpeakwidth, realpeakprom, realpeakheight, fig]=peakfinder_EHedit_10s_withmodeBL(dataforpeakfinder, Fs, ts)
%{
PEAKFINDER
Plots data.
Finds all peaks above threshold. 
Finds median of subthreshold peaks.
Returns all peaks larger than x medians.
Input: 
    mod_dat1: normalized data
    Fs: Sampling frequency
    ts: time vector
Output: 
    peakInd:Index of the peaks in ts
    peakMag: Magnitude of the peaks
%}

minpeakdist= 0.5; %default was 0.2
minpeakheight = 0.5;
tstart_loc=1;
tstop_loc=length(dataforpeakfinder);
t_ind=[tstart_loc:tstop_loc];
colors.green=[0 150 36]/255;
colors.red=[229 28 35]/255;
colors.blue=[0 0 255]/255;

%Find peak threshold: 2 mads over the median of the normalized data
thr=2;
medval=median(dataforpeakfinder); 
smmad=mad(dataforpeakfinder);
thresh=medval+thr*smmad;

%this is EH edit for a moving threshold; number * Fs is number of seconds
%you want the moving window. trying 10s for now?
moving_medval = zeros(1,length(dataforpeakfinder));
moving_smmad = zeros(1, length(dataforpeakfinder));
moving_medval(1:10*Fs) = median(dataforpeakfinder(1:10*Fs));
moving_smmad(1:10*Fs) = mad(dataforpeakfinder(1:10*Fs));
moving_medval(end-(10*Fs):end) = median(dataforpeakfinder(end-(10*Fs):end));
moving_smmad(end-(10*Fs):end) = mad(dataforpeakfinder(end-(10*Fs):end));

moving_mode = zeros(1,length(dataforpeakfinder));
moving_mode(1:10*Fs) = mode(dataforpeakfinder(1:10*Fs));
moving_mode(end-(10*Fs):end) = mode(dataforpeakfinder(end-(10*Fs):end));


for i = (10*Fs)+1:(length(dataforpeakfinder)-(10*Fs))
    moving_medval(i) = median(dataforpeakfinder((i-10*Fs):(i+10*Fs)));
    moving_smmad(i) = mad(dataforpeakfinder((i-10*Fs):(i+10*Fs)));
    moving_mode(i) = mode(dataforpeakfinder((i-10*Fs):(i+10*Fs)));
end

moving_thresh = zeros(1,length(dataforpeakfinder));
for i = 1:length(dataforpeakfinder)
    moving_thresh(i) = moving_medval(i)+thr*moving_smmad(i);
end



%Plot Data
fig = figure;
plot(ts(tstart_loc:tstop_loc), dataforpeakfinder(tstart_loc:tstop_loc))
hold on
%plot(ts([1 length(normDat_DB_modZ)]),[1 1]*thresh,'color',colors.green);
%EH EDIT for the moving threshold
plot(ts(tstart_loc:tstop_loc), moving_thresh(tstart_loc:tstop_loc),'color',colors.green);


%Find all local maxima (no threshold)
[allpeaks, allpeaklocs, allpeakswidth, allpeaksprominence]=findpeaks(dataforpeakfinder(t_ind),'MinPeakDistance',round(minpeakdist*Fs), 'MinPeakProminence',.5,'Annotate','extents');
allpeaklocs=allpeaklocs+tstart_loc-1;

%ok, so my logic is that we first need to eliminate all the data points
%that are over that moving threshold I just made, then get a new median/mad
%that excludes our outliers, and this is the real moving threshold for the
%peaks. In order to do this, I make a new array that sets all values above
%the moving threshold for that time point as NaNs, then use functions that
%exclude NaNs to come up with a new moving threshold.


sub_movingthresh=zeros(1,length(dataforpeakfinder));
for q=1:length(dataforpeakfinder)
    if dataforpeakfinder(q)>moving_thresh(q)
        sub_movingthresh(q)=NaN;
    else sub_movingthresh(q)=dataforpeakfinder(q);
    end
end
sub_movingthresh=transpose(sub_movingthresh);

%Getting official peak threshold: 3 MADs above median, excluding values
%above the first threshold
newthr=3.33;
newmedval=nanmedian(sub_movingthresh); 
newsmmad=mad(sub_movingthresh);
newthresh=newmedval+newthr*newsmmad;

new_moving_medval = zeros(1,length(sub_movingthresh));
new_moving_smmad = zeros(1, length(sub_movingthresh));
new_moving_medval(1:10*Fs) = nanmedian(sub_movingthresh(1:10*Fs));
new_moving_smmad(1:10*Fs) = mad(sub_movingthresh(1:10*Fs));
new_moving_medval(end-(10*Fs):end) = nanmedian(sub_movingthresh(end-(10*Fs):end));
new_moving_smmad(end-(10*Fs):end) = mad(sub_movingthresh(end-(10*Fs):end));

for j = (10*Fs)+1:(length(sub_movingthresh)-(10*Fs))
    new_moving_medval(j) = nanmedian(sub_movingthresh((j-10*Fs):(j+10*Fs)));
    new_moving_smmad(j) = mad(sub_movingthresh((j-10*Fs):(j+10*Fs)));
end

new_moving_thresh = zeros(1,length(sub_movingthresh));
for j = 1:length(sub_movingthresh)
    new_moving_thresh(j) = new_moving_medval(j)+newthr*new_moving_smmad(j);
end

realpeakMag=[];
realpeakInd=[];
realpeakwidth=[];
realpeakprom=[];
realpeakheight=[];

for r=1:length(allpeaks)
    if allpeaks(r)>=new_moving_thresh(allpeaklocs(r))
        realpeakMag(length(realpeakMag)+1)=allpeaks(r);
        realpeakInd(length(realpeakInd)+1)=allpeaklocs(r);
        realpeakwidth(length(realpeakwidth)+1)=allpeakswidth(r);
        realpeakprom(length(realpeakprom)+1)=allpeaksprominence(r);
        realpeakheight(length(realpeakheight)+1)=(allpeaks(r)-new_moving_medval(allpeaklocs(r)));
    end
end

%peakLoc=ts(peakInd);
realpeakloc=ts(realpeakInd);

%Plot peaks
plot(ts(realpeakInd),realpeakMag,'o','markersize',5,'linewidth',2,'color',colors.red);
plot(ts(tstart_loc:tstop_loc), new_moving_thresh(tstart_loc:tstop_loc),'color',colors.red);
plot(ts(tstart_loc:tstop_loc), new_moving_medval(tstart_loc:tstop_loc),'color',colors.blue);

hold off

end