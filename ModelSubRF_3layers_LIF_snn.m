
% updated by JSS
% 2021.4.25  6 - 2- 1 LIF

% a model for subRF 
% there are some dummy variables never used by the code
% 2014.04.23

function ModelSubRF_3layers_LIF_snn(name1,Weight,Weight2,SEED,V_reset,V_e,V_th)

global RefreshRate; %Stimulus refresh rate (Stim frames per second)
global pStim;
global pRF;

global CB;
global NL;
global pWorld;
RefreshRate = 25;  %pRF.maxt = 500; % ms; the filter time scale

Rm = 1; % membrane resistance
tau_m = 10; % membrane time constant
nAmp=0.01;
dt = 1; %1ms
rng(1000);
NL = 1;
comp = -1;

%set up 'world' parameters
pWorld.visAngle = 10;     % horizontal extent of visual angle
pWorld.dur = 10^5;
pWorld.rf = RefreshRate;           % Monitor refresh rate (Hz)
pWorld.dt  = 1/pWorld.rf;
pWorld.n = [1,1,pWorld.dur/pWorld.dt+1];   %[ny, nx, nt] in pixels
pWorld.background = 0;     %background color: [-1,1]
pWorld = initWorld(pWorld); %get the x and y meshrid matrices

% Linear RF in Space and Time
clear pRF
%temporal parameters 
dummy = ceil(1000/RefreshRate); % ms; the bin size to count spikes
pRF.maxt = 800; % ms; the filter time scale
pRF.nt = pRF.maxt/dummy;
pRF.x = fliplr(abs(-pRF.maxt+dummy:dummy:0))/1000;
pRF.xi = 0:pRF.x(2)/5:max(pRF.x);
pRF.t = -fliplr(pRF.xi);

pRF.n = [4,5]; %n cascades
pRF.theta = [.05,.05]; %sec
pRF.a = [.05,.05];
pRF.delay = 2*pWorld.dt; % frames
RF = make1DOG(pWorld,pRF);
tmSTA = squeeze(RF);

pStim.NumImages = 9*10^5;
pStim.xPix = 8;
pStim.yPix = 8;
pStim.Nx = pStim.xPix;
pStim.Ny = pStim.yPix;
pStim.imSize = 4;

% spSTA
spSTA1 = zeros(pStim.xPix, pStim.yPix);
spSTA2 = zeros(pStim.xPix, pStim.yPix);
spSTA3 = zeros(pStim.xPix, pStim.yPix);
spSTA4 = zeros(pStim.xPix, pStim.yPix);
spSTA5 = zeros(pStim.xPix, pStim.yPix);
spSTA6 = zeros(pStim.xPix, pStim.yPix);
spSTA1(2+[1:pStim.imSize/2], 1+[1:pStim.imSize/2]) = 1;
spSTA2(2+[1:pStim.imSize/2], 1+pStim.imSize/2+[1:pStim.imSize/2]) = 1;
spSTA3(2+[1:pStim.imSize/2], 1+pStim.imSize+[1:pStim.imSize/2]) = 1;
spSTA4(2+pStim.imSize/2+[1:pStim.imSize/2], 1+[1:pStim.imSize/2]) = 1;
spSTA5(2+pStim.imSize/2+[1:pStim.imSize/2], 1+pStim.imSize/2+[1:pStim.imSize/2]) = 1;
spSTA6(2+pStim.imSize/2+[1:pStim.imSize/2], 1+pStim.imSize+[1:pStim.imSize/2]) = 1;

spM=[];
spM(1,:) = reshape(spSTA1,1,[]);
spM(2,:) = reshape(spSTA2,1,[]);
spM(3,:) = reshape(spSTA3,1,[]);
spM(4,:) = reshape(spSTA4,1,[]);
spM(5,:) = reshape(spSTA5,1,[]);
spM(6,:) = reshape(spSTA6,1,[]);
%for checkboard
pStim.blackwhite = false;
pStim.meanintensity = 0.5;
pStim.contrast = 0.5;
pStim.seed = SEED;

%Binary checkerboard
CB = ran1(pStim.seed, pStim.Nx*pStim.Ny*pStim.NumImages);
CB = reshape(CB,pStim.Nx*pStim.Ny,pStim.NumImages);
CB(CB>=0.5) = 1;
CB(CB<0.5)  = 0;
CB = CB - 0.5;

Nsub = size(spM,1); 
sum_RF = zeros(Nsub,pStim.NumImages);
tmM=zeros(Nsub,20);
for tem=1:Nsub
    tmM(tem,:)=tmSTA;
end
spSTA = [];
for i=1:Nsub
	temp = reshape(spM(i,:),[],1);
	temp = repmat(temp,1,pStim.NumImages);
	stim = CB .* temp; 
	stim = sum(stim,1);
	a = sameconv(stim',tmM(i,:)');
	% nonlinearity here
    if NL == 1
        a = NL1(a);
    elseif NL == 2
        a = NL2(a,comp);
    elseif NL == 3
        a = NL3(a,comp);
    end
    np = 2;   % should be small, like 2, can not be large
    abig = repmat(a'*pWorld.dt,np,1); % make Poisson spike train
    subsp = ceil(sum(rand(size(abig))<abig,1)');
    sub_SP(i,:) = subsp';
end

gclist=[1,4,2,5;2,5,3,6];
for ig=1:2
    sub_SPt= sub_SP(gclist(ig,:),:);
    output = Weight * sub_SPt; 
    frameT=1000/RefreshRate;
    insp=zeros(1,1000/RefreshRate*length(output));
    insp(1000/RefreshRate/2:1000/RefreshRate:end)=output;
    T = 0:1:length(insp); 
    Vm(1) = V_e;
    Vm(2) = V_e;
    Im=zeros(1,length(T));
    tau=10;
    fdt=1-(dt/tau);
    spk = num2cell(ones(1));    % cell array of spikes
    
    for t=3:length(T)
        Im(t) = (2*fdt*Im(t-1)) - ((fdt^2)*Im(t-2)) + (((dt/tau)^2) * insp(t-2));
      	% integration
     	if Vm(t-1) >= V_th
          	Vm(t) = V_reset;
        else
            Vm(t) = Vm(t-1) + dt * ( -(Vm(t-1) - V_e) + Im(t) * Rm +nAmp.*randn(1) )/tau_m ;
        end
       	% spiking
      	isSpk = Vm(:,t)>=V_th;
       	if any(isSpk)
          	% append spike times to spk
           	spk(isSpk) = cellfun( @(u,v) [u v],spk(isSpk),num2cell(isSpk(isSpk)'*t),'unif',false );
           	% reset membrane voltage
         	Vm(isSpk,t) = V_reset;
       	end   
    end

    Vm=Vm(2:end);
    spkTime=spk{1};
    spkTime=spkTime(2:end);
    frameT=1000/RefreshRate;
    ftimes=1:frameT:length(Vm);
    spike = histc( spkTime,ftimes )';  % 25Hz
    spike(1:pRF.nt)=0;
    sub_SP2(ig,:)=spike;
end 

    output2 = Weight2 * sub_SP2; 
    insp=zeros(1,1000/RefreshRate*length(output2));
    insp(1000/RefreshRate/2:1000/RefreshRate:end)=output2;
    T = 0:1:length(insp); 
    Vm(1) = V_e;
    Vm(2) = V_e;
    Im=zeros(1,length(T));
    tau=10;
    fdt=1-(dt/tau);
    spk = num2cell(ones(1));    % cell array of spikes
    
    for t=3:length(T)
        Im(t) = (2*fdt*Im(t-1)) - ((fdt^2)*Im(t-2)) + (((dt/tau)^2) * insp(t-2));
      	% integration
     	if Vm(t-1) >= V_th
          	Vm(t) = V_reset;
        else
            Vm(t) = Vm(t-1) + dt * ( -(Vm(t-1) - V_e) + Im(t) * Rm +nAmp.*randn(1) )/tau_m ;
        end
       	% spiking
      	isSpk = Vm(:,t)>=V_th;
       	if any(isSpk)
          	% append spike times to spk
           	spk(isSpk) = cellfun( @(u,v) [u v],spk(isSpk),num2cell(isSpk(isSpk)'*t),'unif',false );
           	% reset membrane voltage
         	Vm(isSpk,t) = V_reset;
       	end   
    end;

    Vm=Vm(2:end);
    spkTime=spk{1};
    spkTime=spkTime(2:end);
    ftimes=1:frameT:length(Vm);
    spike = histc( spkTime,ftimes )';  % 25Hz
    spike(1:pRF.nt)=0;

    spklist = spike(spike>0);

    [sta, ~, SS]= simpleSTC(CB',spike,pRF.nt);
    sta = sta./norm(sta);
    sta = reshape(sta,pRF.nt,pStim.Nx, pStim.Nx);
    nt = pRF.nt;
    nx = pStim.Nx;
    ny = pStim.Ny;
    
    % singular value decomposition
    temp=[];
    [temp,sigma,sp] = svd(reshape(sta,size(sta,1),[]),'econ');
    sigma = diag(sigma);
    svdsta = cell(size(sigma));
    for m=1:length(svdsta)
        svdsta{m} = reshape(temp(:,m)*sigma(m)*sp(:,m)',size(sta));
    end
    tmSVD = temp(:,1);
    if abs(min(tmSVD))==max(abs(tmSVD))
        tmSVD=-tmSVD;
        labelF=-1;
    else
        labelF=1;
    end
    for n=1:1
        % spatial SVD component
        spSVD = reshape(sp(:,n)',size(svdsta{n},2),size(svdsta{n},3));
        spSVD=labelF*spSVD;
        % visualize spatial SVD component
        % data range
        datarange = [-0.5 0.5];
        if ~diff(datarange), datarange = mean(datarange)+[-0.1,0.1]; end
    end
    
    save([name1,'Data.mat'],'CB','pStim','spkTime','spike','SS','sub_SP','sub_SP2','spklist','nt','nx','ny','tmSVD','spSVD','spM','tmM','sta','comp','-v7.3');
end


function [x] = NL1(x)
x(x<0) = 0;
end

%sublinear
function [x] = NL2(x,comp)
fun1 = @(x) 0.5 - sqrt(0.5.^2-x.^2);
ind = find(x<0);
x = fun1(x);
x(ind) = -x(ind);

if comp<0
    x(x>0) = 0;
elseif comp>0
    x(x<0) = 0;
end
end

%superlinear
function [x] = NL3(x,comp)

fun1 = @(x) sqrt( 0.5.^2 - (x-0.5).^2 );
fun2 = @(x) -sqrt( 0.5.^2 - (x+0.5).^2 );
if comp<0
    x(x>0) = 0;
    x = fun2(x);
elseif comp>0
    x(x<0) = 0;
    x = fun1(x);
end
end

