function RippleysZerosArea
%This code created to calculate and plot Rippley's K,L,H functions
%Version 1.0 11:55 August 9 2017: General functionality
%Version 2.0 13:32 September 7 2017: Calculate zeros, peaks, area
%Version 2.1 12:07 September 11 2017: Repair dimensional integration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all
clc

% prompt = 'What cellular disances would you like to evauluate over?(enter a matrix)\n';
% distpick = input(prompt);
% prompt = 'What frame(s) would you like to evaluate over?\n';
% frameset = input(prompt);
% loadpath = uigetdir('Z:\ENG_BBCancer_Shared','Where would you like to load the data from?');
% prompt = 'What well(s) would you like to evalaute?\n';
% wells = input(prompt);
% prompt = 'Would you like to save the image? (1 for yes)\n';
% saveind = input(prompt);
% if saveind == 1
%     savefolder = uigetdir('Z:\ENG_BBCancer_Shared','Where would you like to save the images?');
% end
distpick = [1:400];
distinc = 1;
frameset = [1:20:241];
loadpath = 'Z:\ENG_BBCancer_Shared\group\0Zach\Cluster data\EGF (E6) Data';
% wells = [2 3 4 5 6 19 20 21 22 23 24 49 50 51 52 53 54 67 68 69 70 71 72];
wells = [2 3 4 5];
saveind = 0;
areacell = cell(1,max(wells));
wellcounter = 1;

for pickwell = wells
%set up color matp
cmap = cbrewer('seq','YlOrRd',numel(frameset)+1);
colormap(cmap);
figure;
% filename = strcat('w',num2str(pickwell),'-SCALEDalltrackfile.mat');
filename = strcat(loadpath,'\egf(E6)newwell',num2str(pickwell),'.mat');
% load super track file
load(filename)

%preallocate kmat matrix
lmat = zeros(numel(distpick),numel(frameset));
hmat = zeros(numel(distpick),numel(frameset));
maxmat = zeros(numel(frameset),2);
zeromat = zeros(numel(frameset),1);
area = zeros(numel(frameset),2);

%loop through the respective time frames
framecounter = 0;
for frame = frameset
    framecounter = framecounter + 1;
    %pull out the values that are recoded
    Vals = find(~isnan(storeX(:,frame))); %#ok<NODEF>
    %adjust x and y accordingly
    tempX = storeX(Vals,frame);
    tempY = storeY(Vals,frame);
    %preallocate K matrix

    K = zeros(numel(distpick),1);
    L = zeros(numel(distpick),1);
    H = zeros(numel(distpick),1);
    plotdist = distpick;
    b = 1;
    for distance = distpick
        K(b) = RipleysK([tempX tempY],distance,[0 1000 0 1000],0);
        L(b) = sqrt(K(b)/pi);
        H(b) = L(b) - distance;
        b = b+1;
    end
    %find maximum value of the H function
    [~,Ixy] = max(H);
    maxmat(framecounter,:) = [plotdist(Ixy),H(Ixy)];
    
    %locate first negative value to the left of the local maxima
    checkmat = fliplr(0:maxmat(framecounter,1));
    for uu = checkmat
        if H(uu) < 0
            tempzero = uu+1;
            break
        end
    end
    
    zeromat(framecounter) = tempzero;
    hrest = H(tempzero:maxmat(framecounter,1));
    area(framecounter,:) = [frame trapz(hrest,distinc)];

    lmat(:,framecounter)=L;
    hmat(:,framecounter)=H;
    ax1 = subplot(2,1,1);
    %plot the function
%     ylim([0 300]);
    plot(plotdist,L,'-','Linewidth',2,'Color',cmap(framecounter,:));
    hold on
    pbaspect(ax1,[3 2 2]);
    ax2 = subplot(2,1,2);
%     ylim([-100 100]);
    plot(plotdist,H,'-','Linewidth',2,'Color',cmap(framecounter,:));
    pbaspect(ax2,[3 2 2]);
    hold on
    fprintf('frame %d is now complete\n',frame)
end
    areacell{1,pickwell} = area;
    wellcounter = wellcounter + 1;

%plotting
subplot(2,1,1);

plot(plotdist,plotdist,'color',[0 0 1],'linewidth',1.5);
xlabel('Distance in Microns');
ylabel('Ripleys L Function');
titlename = strcat('Well',num2str(pickwell),'Plot');
title(titlename);

subplot(2,1,2);
ylim([-100 100])
plot(plotdist,plotdist.*0,'color',[ 0 0 1],'linewidth',1.5);
xlabel('Distance in Microns');
ylabel('Ripleys H Function');
if saveind == 1
    savefilename = strcat(savefolder,'\W',num2str(pickwell),'ripleysfunctionplot2000');
    print(gcf,savefilename,'-dpng','-r2000')
    print(gcf,savefilename,'-depsc','-r2000')
end
end
beep
% save('IntegralDataRippleys91117','areacell')
excelmat = zeros(numel(frameset)+1,numel(wells)*2);
excelmatind = 1;

%compile data in one spreadsheet
for ii = 1:numel(areacell)
   if numel(areacell{ii}) == 0
   else
       excelmat(1,excelmatind) = ii;
       excelmat(2:end,excelmatind:excelmatind+1) = areacell{ii};
       excelmatind = excelmatind + 2;
   end
end



end

function K = RipleysK(locs,dist,box,method)
% RipleysK: Calculate K statistic
% 
% K = RipleysK(locs,dist, box,method) calculates G, the K statistic at each 
% distance for the data with x-y coordinates given by locs, and the
% bounding rectangle given by box=[minX maxX minY maxY].
% If method=0, no edge correction is applied.
% If method=1, points are only used if they are at least h units away from
% the edge.
%
% Note: The L statistic may be calculated from the K statistic as follows: 
%   L = sqrt(K/pi)-h;
%   

if nargin<4, method=1; end
[N,k] = size(locs);
if k~=2, error('locs must have two columns'); end
rbox = min([locs(:,1)'-box(1);box(2)-locs(:,1)';locs(:,2)'-box(3); box(4)-locs(:,2)'] );
% rbox is distance to box

DX = repmat(locs(:,1),1,N)-repmat(locs(:,1)',N,1);
DY = repmat(locs(:,2),1,N)-repmat(locs(:,2)',N,1);
DIST = sqrt(DX.^2+DY.^2);
DIST = sort(DIST);


if method==1
K = zeros(length(dist),1);
for k=1:length(K)
    I = find(rbox>dist(k));
    if length(I)>0 %#ok<ISMT>
        K(k) = sum(sum(DIST(2:end,I)<dist(k)))/length(I);
    end
end
elseif method==0
    K = zeros(length(dist),1);
    for k=1:length(K)

        K(k) = sum(sum(DIST(2:end,:)<dist(k)))/N;
    end
end

lambda = N/((box(2)-box(1))*(box(4)-box(3)));
K = K/lambda;
end