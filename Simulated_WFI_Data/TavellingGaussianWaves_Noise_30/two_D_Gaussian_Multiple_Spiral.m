function two_D_Gaussian_Multiple_Spiral

numberOfGridPoints = 128; % number of pixels in each x and y directions
frameRate = 150; % hz
frameCounter = 1;
durations{1} = [0.1 0.2 0.1];
durations{2} =  [0.1 0.2 0.1];
durations{3} =  [0.1 0.1 0.1];
startTimeOfSources = [0 0.1 0.15];
for ii = 1:length(durations)
    allDurations(ii) = sum(durations{ii});
    tempVar(ii) = allDurations(ii) + startTimeOfSources(ii);
end
% totalTime = max(durations1 + startTimeOfSources(1),allDurations(2) + startTimeOfSources(2));
totalTime = max(tempVar);

% durations1 = sum(durations{1});
% durations2 = sum(durations{2});
% totalTime = max(durations1 + startTimeOfSources(1),allDurations(2) + startTimeOfSources(2));
% totalTime = durations1;
totalFrames =round( totalTime * frameRate);
startFrameOfSources = round(startTimeOfSources * frameRate) + 1;
numberOfSources = 3;

for ii = 1:numberOfSources
    numberOfFrames{ii} = round(durations{ii} * frameRate);
end
for ii = 1:numberOfSources
    totalFramesSources(ii) = sum(numberOfFrames{ii});
end

visualization = [0 0 0];
[gP Va1] = source1(numberOfFrames{1});
frames{1} = generateFrames(numberOfGridPoints,gP,visualization(1));
[gP Va2] = source2(numberOfFrames{2});
frames{2} = generateFrames(numberOfGridPoints,gP,visualization(2));
[gP Va3] = source3(numberOfFrames{3});
frames{3} = generateFrames(numberOfGridPoints,gP,visualization(3));


allFrames1 = zeros([size(frames{1}(:,:,1)) totalFrames]);
allFrames1(:,:,startFrameOfSources(1):(startFrameOfSources(1)+totalFramesSources(1)-1)) = frames{1};
% zP = zeros(1,totalFrames-totalFramesSources(1));
% Va1 = [Va1 zP];

allFrames2 = zeros([size(frames{1}(:,:,1)) totalFrames]);
allFrames2(:,:,startFrameOfSources(2):(startFrameOfSources(2)+totalFramesSources(2)-1)) = frames{2};
% zP = zeros(1,totalFrames-totalFramesSources(2));
% Va2 = [zP Va2];

allFrames3 = zeros([size(frames{1}(:,:,1)) totalFrames]);
allFrames3(:,:,startFrameOfSources(3):(startFrameOfSources(3)+totalFramesSources(3)-1)) = frames{3};

allFrames = allFrames1+allFrames2+allFrames3;
ImgSeq = allFrames;
allFrames = addNoiseToImgSeq(ImgSeq,30);
for ii = 1:totalFrames
    figure(1);clf;
    imagesc(allFrames(:,:,ii));
    axis off;
    axis equal;
    title(sprintf('%d',ii));
end

ImgSeq = allFrames;
fr11 = frames{1}(:,:,16);
fr12 = frames{1}(:,:,45);
fr21 = frames{2}(:,:,16);
fr22 = frames{2}(:,:,45);

figFrameS1 = (fr11 + fr12)/2;
figFrameS2 = (fr21 + fr22)/2;
figFrameS3 = frames{3}(:,:,15);
save('ImgSeq.mat','ImgSeq','Va1','Va2','figFrameS1','figFrameS2','figFrameS3','startFrameOfSources','totalFramesSources');

totalFrames  = size(allFrames,3)


function frameSeq = generateFrames (numberOfGridPoints,gP,v)
allxo = gP(1,:);
allyo = gP(2,:);
allA = gP(3,:);
allsigmaX = gP(4,:);
allsigmaY = gP(5,:);
allRA = gP(6,:);

for ii = 1:length(allxo)
    gaussParam.center = [allxo(ii) allyo(ii)];
    gaussParam.amplitude = allA(ii);
    gaussParam.sigma = [allsigmaX(ii) allsigmaY(ii)];
    frameSeq(:,:,ii) = generateRotatedFrame(numberOfGridPoints,gaussParam,allRA(ii));
    if v
        figure(1);clf;
        imagesc(frameSeq(:,:,ii));
        axis equal
        axis off
    end
end

function frame = generateRotatedFrame (numberOfGridPoints,gaussParam,gRA)
% this function generates a 2D gaussian provided the gauss param and number
% of grid points
x0 = gaussParam.center(1) ; y0 = gaussParam.center(2) ;
A = gaussParam.amplitude;
sigma_x = gaussParam.sigma(1); sigma_y = gaussParam.sigma(2);
temp = linspace(ceil(-numberOfGridPoints/2)+1,ceil(numberOfGridPoints/2),numberOfGridPoints);
[X,Y] = meshgrid(temp,temp);
a = cos(gRA)^2/(2*sigma_x^2) + sin(gRA)^2/(2*sigma_y^2);
b = -sin(2*gRA)/(4*sigma_x^2) + sin(2*gRA)/(4*sigma_y^2);
c = sin(gRA)^2/(2*sigma_x^2) + cos(gRA)^2/(2*sigma_y^2);

frame = A*exp( - (a*(X-x0).^2 - 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
frame  = flipud(frame);


function [val Va] = source1(numberOfFrames)

sigmaSource = [15];
sM =[1] ;
allxo = -20 * ones(1,numberOfFrames(1)); 
allyo =-30*ones(1,numberOfFrames(1));
allA = 1*ones(1,numberOfFrames(1));
allsigmaX = linspace(0.01,sigmaSource(1),numberOfFrames(1)); 
allsigmaY = sM(1)*allsigmaX;
allRA = zeros(1,numberOfFrames(1));
Va = zeros(1,numberOfFrames(1));
xo = allxo(end);
yo = allyo(end);
radius = sqrt(xo^2 + yo^2);
startingAngle = atan2(yo,xo);%-pi/2 + (-pi/4);
vf = 0.5/numberOfFrames(2);
tFrames = linspace(0,numberOfFrames(2),numberOfFrames(2));
w = sin(2*pi*vf*tFrames);
theta = -cos(2*pi*vf*tFrames)+startingAngle + 1;
% theta = (1/(2*pi*vf))*(-cos(2*pi*vf*tFrames)+startingAngle*2*pi*vf + 1);
Va = [Va 0 radius * (diff(theta)./diff(tFrames))];
angles = theta;
A = 1;
RA = linspace(0,0,numberOfFrames(2));
for ii = 1:length(angles)
    thisAngle = angles(ii);
    xc = round(radius * cos(thisAngle));
    yc = round(radius * sin(thisAngle));
    allxo = [allxo xc];
    allyo = [allyo yc];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaSource(1)];
    allsigmaY = [allsigmaY sM(1)*sigmaSource(1)];
    allRA = [allRA RA(ii)];
end
xo = allxo(end); yo = allyo(end);  % center of gaussian
A = 1;
sigmaX = linspace(sigmaSource(1),0.01,numberOfFrames(3));
sigmaY = sM(1)*sigmaX;
for ii = 1:length(sigmaX)
    allxo = [allxo xc];
    allyo = [allyo yc];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaX(ii)];
    allsigmaY = [allsigmaY sigmaY(ii)];
    allRA = [allRA RA(end)];
end
tVa = zeros(1,numberOfFrames(3));
Va = [Va tVa];
val = [allxo;allyo;allA;allsigmaX;allsigmaY;allRA];

function [val Va] = source2(numberOfFrames)
sigmaSource = 6;
sM = 2;
allxo = 30 * ones(1,numberOfFrames(1)); 
allyo =30*ones(1,numberOfFrames(1));
allA = 1*ones(1,numberOfFrames(1));
allsigmaX = linspace(0.01,sigmaSource(1),numberOfFrames(1)); 
allsigmaY = sM(1)*allsigmaX;
allRA = zeros(1,numberOfFrames(1));
Va = zeros(1,numberOfFrames(1));
xo = allxo(end);
yo = allyo(end);
radius = sqrt(xo^2 + yo^2);
startingAngle = atan2(yo,xo);%-pi/2 + (-pi/4);
vf = 0.5/numberOfFrames(2);
tFrames = linspace(0,numberOfFrames(2),numberOfFrames(2));
w = sin(2*pi*vf*tFrames);
theta = -cos(2*pi*vf*tFrames)+startingAngle + 1;
Va = [Va 0 radius * (diff(theta)./diff(tFrames))];
angles = theta;
A = 1;
RA = linspace(0,pi/6,numberOfFrames(2));
for ii = 1:length(angles)
    thisAngle = angles(ii);
    xc = round(radius * cos(thisAngle));
    yc = round(radius * sin(thisAngle));
    allxo = [allxo xc];
    allyo = [allyo yc];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaSource(1)];
    allsigmaY = [allsigmaY sM(1)*sigmaSource(1)];
    allRA = [allRA RA(ii)];
end
xo = allxo(end); yo = allyo(end);  % center of gaussian
A = 1;
sigmaX = linspace(sigmaSource(1),0.01,numberOfFrames(3));
sigmaY = sM(1)*sigmaX;
for ii = 1:length(sigmaX)
    allxo = [allxo xc];
    allyo = [allyo yc];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaX(ii)];
    allsigmaY = [allsigmaY sigmaY(ii)];
    allRA = [allRA RA(end)];
end
val = [allxo;allyo;allA;allsigmaX;allsigmaY;allRA];
tVa = zeros(1,numberOfFrames(3));
Va = [Va tVa];


function [val Va] = source3(numberOfFrames)

sigmaSource = [10];
sM =[1] ;
allxo = -40 * ones(1,numberOfFrames(1)); 
allyo =-10*ones(1,numberOfFrames(1));
allA = 1*ones(1,numberOfFrames(1));
allsigmaX = linspace(0.01,sigmaSource(1),numberOfFrames(1)); 
allsigmaY = sM(1)*allsigmaX;
allRA = zeros(1,numberOfFrames(1));
tempxo = linspace(-40,-20,numberOfFrames(2));
yo = allyo(end);
A = allA(end);
sigmaX = allsigmaX(end);
sigmaY = allsigmaY(end);
for ii = 1:length(tempxo)
    allxo = [allxo tempxo(ii)];
    allyo = [allyo yo];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaX];
    allsigmaY = [allsigmaY sigmaY];
    allRA = [allRA 0];
end
RA = 0;
xo = allxo(end); yo = allyo(end);  % center of gaussian
A = 1;
sigmaX = linspace(sigmaSource(1),0.01,numberOfFrames(2));
sigmaY = sM(1)*sigmaX;
for ii = 1:length(sigmaX)
    allxo = [allxo xo];
    allyo = [allyo yo];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaX(ii)];
    allsigmaY = [allsigmaY sigmaY(ii)];
    allRA = [allRA RA(end)];
end
val = [allxo;allyo;allA;allsigmaX;allsigmaY;allRA];
Va = 0;