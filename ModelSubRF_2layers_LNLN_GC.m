% a model for subRF 

function ModelSubRF_2layers_LNLN_GC(pathname,Weight)

    global RefreshRate; %Stimulus refresh rate (Stim frames per second)
    global pStim;
    global pRF;
    RefreshRate = 30;

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
    %tmSTA 1-4 different amplitudes  same delay  
    %tmSTA 5-8 same amplitudes different delay  and tmSTA7=tmSTA3;
    dummy = ceil(1000/RefreshRate); % ms; the bin size to count spikes
    pRF.maxt = 680; % ms; the filter time scale
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


    %pStim.NumImages = 12*10^6;
    pStim.NumImages = 10^5;
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
    %spSTA(2+[1:pStim.imSize], 2+[1:pStim.imSize]) = 1;
    spSTA1(2+[1:pStim.imSize/2], 2+[1:pStim.imSize/2]) = 1;
    spSTA2(2+[1:pStim.imSize/2], 2+pStim.imSize/2+[1:pStim.imSize/2]) = 1;
    spSTA3(2+pStim.imSize/2+[1:pStim.imSize/2], 2+[1:pStim.imSize/2]) = 1;
    spSTA4(2+pStim.imSize/2+[1:pStim.imSize/2], 2+pStim.imSize/2+[1:pStim.imSize/2]) = 1;

    spM=[];
    spM(1,:) = reshape(spSTA1,1,[]);
    spM(2,:) = reshape(spSTA2,1,[]);
    spM(3,:) = reshape(spSTA3,1,[]);
    spM(4,:) = reshape(spSTA4,1,[]);

    %for checkboard
    pStim.blackwhite = false;
    pStim.meanintensity = 0.5;
    pStim.contrast = 0.5;
    pStim.seed = -1000;

    %Bineary checkerboard
    CB = ran1(pStim.seed, pStim.Nx*pStim.Ny*pStim.NumImages);
    CB = reshape(CB,pStim.Nx*pStim.Ny,pStim.NumImages);
    CB(CB>=0.5) = 1;
    CB(CB<0.5)  = 0;
    CB = CB - 0.5;

    Nsub = 4;
    sum_RF = zeros(Nsub,pStim.NumImages);
    tmM=zeros(4,20);
    tmM(1,:)=tmSTA;
    tmM(2,:)=tmSTA;
    tmM(3,:)=tmSTA;
    tmM(4,:)=tmSTA;
    spSTA = [];
    for i=1:Nsub
        if i==1
            spSTA = spSTA1;
        elseif i==2
            spSTA = spSTA2;
        elseif i==3
            spSTA = spSTA3;
        elseif i==4
            spSTA = spSTA4;
        end
        temp = reshape(spSTA,[],1);
        temp = repmat(temp,1,pStim.NumImages);
        stim = CB .* temp; 
        stim = sum(stim,1);
        a = sameconv(stim',tmSTA);
        % nonlinearity here
        if NL == 1
            a = NL1(a);
        elseif NL == 2
            a = NL2(a,comp);
        elseif NL == 3
            a = NL3(a,comp);
        end
        sum_RF(i,:) = a;
    end

    %get data
    output = Weight * sum_RF; 
    pModel.a = 2;  %scale fac tor
    pModel.b = 10;         %exponent
    pModel.Vrest =2.5;     %threshold
    r = 1*output-pModel.Vrest;
    r(r<0) = 0;
    r = 2*r;
    np = 2;   % should be small, like 2, can not be large
    rbig = repmat(r*pWorld.dt,np,1); % make Poisson spike train
    sp = ceil(sum(rand(size(rbig))<rbig,1)');
    sp2 = sp/(np*pWorld.dt);
    spk = sp;
    sp(1:pRF.nt) = 0;
    spike=sp;
    spklist = sp(sp>0);
    nt = pRF.nt;
    nx = pStim.Nx;
    ny = pStim.Ny;
    
    save([pathname,'GCData_model.mat'],'CB','pStim','spike','nt','nx','ny','spM','tmM','comp','-v7.3');
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


function [output] = modelSpk(p, pdata)
data = pdata{1};
mean = pdata{2};
sigfun = pdata{3};
dummy = repmat(p',1,size(data,2))';
X = sum(dummy.*mean,2);
[Vm_IM, sp_IM] = NL_model(X, data);
x = Vm_IM;
y = sp_IM;
xx = diff(x);
xx = sort([x x(1:end-1)+xx/2 ]);
yy = spline(x, y, xx);
pp = mineFit(sigfun,xx,yy,[1 2 0.5]);
output = sigfun(pp, X);
end

