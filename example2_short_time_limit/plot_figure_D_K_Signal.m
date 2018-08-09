

root = '/Users/magda/Documents/GitHub/monte-carlo-simulation-recipes/example2_short_time_limit';

%% Have a look of the microstructure
% The lookup table A saves two axon labels in one integer. If the first
% and the second axon labels are ax1 and ax2, ax1 = mod(A,Nmax), and ax2 =
% floor(A/Nmax).
% Other parameters:
%   voxelSize: voxel size of the look up table A in µm
%   Nax: # axons
%   rCir: axon's inner radius
%   gratio: g-ratio, the ratio of inner to outer radius
%   [xCir,yCir]: axon's center position

root_packing = fullfile(root,'packing');

A = load(fullfile(root_packing,'phantom_APix.txt')); 
Nmax = load(fullfile(root_packing,'phantom_Nmax.txt'));
voxelSize = load(fullfile(root_packing,'phantom_res.txt'));

Nax = load(fullfile(root_packing,'phantom_NAx.txt'));
gratio = load(fullfile(root_packing,'phantom_gratio.txt')); 
rCir = load(fullfile(root_packing,'phantom_rCir.txt'));
xCir = load(fullfile(root_packing,'phantom_xCir.txt'));
yCir = load(fullfile(root_packing,'phantom_yCir.txt'));

% Plot the microstructure and the lookup table
figure;
subplot(121);
for i = -1:1
    for j = -1:1
        viscircles([xCir(:)+i,yCir(:)+j],rCir(:));
    end
end
xlim([0 1]); ylim([0 1])
pbaspect([1 1 1]); box on
title('Axon Packing','interpreter','latex','fontsize',20)
set(gca,'xtick',[],'ytick',[])

subplot(122);
imagesc(rot90(A),[0,Nax]);
pbaspect([1 1 1]); axis off
title('Lookup table','interpreter','latex','fontsize',20)

set(gcf,'unit','inch','position',[0 0 12 5])

%% simulation parameter

clear tline C
fileID = fopen(fullfile(root,'sim_para'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

dt=str2num(C{1});       % time step
TN_max=str2num(C{2});   % total time step
Np=str2num(C{3});       % # of particles
TN_num=str2num(C{4});   % # of recorded time step
Din=str2num(C{5});      % diffusivity in IAS
Dout=str2num(C{6});     % diffusivity in EAS
kappa=str2num(C{7});    % permeability of inner and outer myelin sheath
Dmr=str2num(C{8});      % diffusivity within myelin, perpendicular
Dmc=str2num(C{9});      % diffusivity within myelin, circumferential
b=str2num(C{10});       % bvalue

% dx^2

clear tline C
fileID = fopen(fullfile(root,'dx2_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear dx2
for i=1:n-1
    dx2(i,:)=str2num(C{i});
end

%dy^2

clear tline C
fileID = fopen(fullfile(root,'dy2_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear dy2
for i=1:n-1
    dy2(i,:)=str2num(C{i});
end

% dx^4

fileID = fopen(fullfile(root,'dx4_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear dx4
for i=1:n-1
    dx4(i,:)=str2num(C{i});
end

%dy^4

clear tline C
fileID = fopen(fullfile(root,'dy4_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear dy4
for i=1:n-1
    dy4(i,:)=str2num(C{i});
end

% dt=7.5e5*1e-9;
TN=linspace(1,TN_max,TN_num).';
TD = (1:TN)*dt;
DCumX=dx2/2/dt./TN;
DCumY=dy2/2/dt./TN;

KCumX=dx4./dx2.^2-3;
KCumY=dy4./dy2.^2-3;

% plot figure D and K of Cumulant

figure;
subplot(221);
hold on
hx=plot((TN*dt),DCumX,'bo'); axis square
set(hx,'markersize',6,'markerfacecolor',[0 0 1])
hy=plot((TN*dt),DCumY,'ro');
set(hy,'markersize',6)
set(gca,'xscale','log')
% rin=1; d0=Din;
% TD=1e-5:0.001:(max(TN)*dt);
% DTheory=rin^2/4./TD;
% ht=plot(1./(TD),DTheory,'-k');
% set(ht,'linewidth',1);

d0=Dout;

SoV=2*pi*sum(r) / (1-pi*sum(r.^2)) / res;
DMitra=d0*(1-4/6/sqrt(pi)*SoV*sqrt(d0*TD));
htm=plot((TD),DMitra,'-.');
set(htm,'linewidth',1,'color',[0 0.6 0]);

% N=20;
% dJ0 = @(x) -besselj(1,x);
% beta0k = @(k) fzero(dJ0, [(k-1)*pi, k*pi]); 
% beta0=zeros(1,N); for k=1:N, beta0(k) = beta0k(k+1); end;
% % beta0=[3.8317    7.0156   10.1735   13.3237   16.4706   19.6159   22.7601   25.9037   29.0468   32.1897];
% beta0 = [NaN beta0(1:N-1)];
% 
% % N=10;
% dJ1 = @(x) besselj(0,x) - besselj(1,x)./(x + 1e-15);
% beta1k = @(k) fzero(dJ1, [(k-1)*pi, k*pi]); 
% beta1=zeros(1,N); for k=1:N, beta1(k) = beta1k(k); end;
% % beta1 = [1.8412    5.3314    8.5363   11.7060   14.8636   18.0155   21.1644   24.3113   27.4571   30.6019];
% 
% % N=10;
% dJ2 = @(x) besselj(1,x) - besselj(2,x)*2./(x + 1e-5);
% beta2k = @(k) fzero(dJ2,[floor((k-1)*pi), k*pi]);
% beta2=zeros(1,N); for k=1:N, beta2(k) = beta2k(k+1); end;
% % beta2=[3.0542    6.7061    9.9695   13.1704   16.3475   19.5129   22.6716   25.8260   28.9777   32.1273];
% beta2 = [NaN beta2(1:N-1)];
% 
% 
% DTheoryExact=rin^2/4./TD;
% for i=1:numel(beta1)
%     DTheoryExact=DTheoryExact-2*rin^2./TD.*exp(-beta1(i)^2*d0.*TD/rin^2)/beta1(i)^2./(beta1(i)^2-1);
% end

fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0,0],...
    'Upper',[Inf,Inf,Inf],...
    'StartPoint',[1 1 1]);
ft = fittype('Di+A*(log(t/tc)+3/2)/(t-3.5e-5/3)',...
    'dependent','Dt','independent','t','coefficients',{'Di','A','tc'},...
    'options',fo);
[curve, gof] = fit(TN(50:end)*dt,DCumX(50:end),ft);

DTheoryExact=curve.Di+curve.A*(log(TD/curve.tc)+3/2)./(TD-dt/3);
hte=plot((TD),DTheoryExact,'--k');
set(hte,'linewidth',1);

% numer=5/64/3 + exp(-beta1(1)^2*d0*TD/rin^2)/beta1(1)^2/(beta1(1)^2-1)*(4/beta1(1)^2-1.5);
% denom=1/4 - 2*exp(-beta1(1)^2*d0*TD/rin^2)/beta1(1)^2/(beta1(1)^2-1);
% for i=2:numel(beta0)
%     numer=numer + exp(-beta0(i)^2*d0*TD/rin^2)/beta0(i)^4 + exp(-beta1(i)^2*d0*TD/rin^2)/beta1(i)^2/(beta1(i)^2-1)*(4/beta1(i)^2-3/2) +...
%         exp(-beta2(i)^2*d0*TD/rin^2)/beta2(i)^2/(beta2(i)^2-4)/2;
%     denom=denom - 2*exp(-beta1(i)^2*d0*TD/rin^2)/beta1(i)^2/(beta1(i)^2-1);
% end
% numer=6*numer;
% denom=(denom).^2;
% KTheoryExact=numer./denom-3;
% KMitra=1e3+0*TD;%numer./(d0*TD/rin^2).^2-3;

KTheoryExact=6*curve.A/curve.Di*(log(TD/curve.tc)+3/2)./(TD-dt/3);

legend([hx hy hte htm],{'$D_x=\langle x^2 \rangle /2t$','$D_y=\langle y^2 \rangle /2t$',...
    '$D_{\rm theory\,exact}$','$D_{\rm Mitra}$'},'interpreter','latex','fontsize',20,'location','southwest')

xlim([0 (max(TN)*dt)]); 
ylim([0 2.1]); box on;
xlabel('${\Delta}$ (ms)','interpreter','latex','fontsize',20);
ylabel('$D_\perp$ ($\mu$m$^2$/ms)','interpreter','latex','fontsize',20);

subplot(222);
hold on
hx=plot(sqrt(TN*dt),KCumX,'bo'); axis square
set(hx,'markersize',6,'markerfacecolor',[0 0 1])
hy=plot(sqrt(TN*dt),KCumY,'ro');
set(hy,'markersize',6)
% ht=refline(0,-0.5);
% set(ht,'color',[0 0 0])
hte=plot(sqrt(TD),KTheoryExact,'--k');
set(hte,'linewidth',1);
% htm=plot(sqrt(TD),KMitra,'-.');
% set(htm,'linewidth',1,'color',[0 0.6 0])

legend([hx hy hte],{'$K_x=\frac{\langle x^4 \rangle}{\langle x^2 \rangle ^2} -3$','$K_y=\frac{\langle y^4 \rangle}{\langle y^2 \rangle ^2} -3$',...
    '$K_{\rm theory\,exact}$'},...
    'location','southeast','interpreter','latex','fontsize',20)

xlim([0 sqrt(max(TN)*dt)]); 
ylim([-2 3]); 
box on;
xlabel('$\sqrt{\Delta}$ (ms)','interpreter','latex','fontsize',20);
ylabel('$K_\perp$','interpreter','latex','fontsize',20);

%

% Sx

clear tline C
fileID = fopen(fullfile(root,'gxSigRe_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear SxR
for i=1:n-1
    SxR(i,:)=str2num(C{i});
end

clear tline C
fileID = fopen(fullfile(root,'gxSigIm_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear SxI
for i=1:n-1
    SxI(i,:)=str2num(C{i});
end

% Sy

clear tline C
fileID = fopen(fullfile(root,'gySigRe_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear SyR
for i=1:n-1
    SyR(i,:)=str2num(C{i});
end

clear tline C
fileID = fopen(fullfile(root,'gySigIm_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear SyI
for i=1:n-1
    SyI(i,:)=str2num(C{i});
end

% DKI

% b=0.1:0.1:10;

clear SigX DSigX KSigX DSigY KSigY
SigX=sqrt(SxR.^2+SxI.^2);
SigY=sqrt(SyR.^2+SyI.^2);

bMax=1;
bIdx=find(b<=bMax);
bLow=b(bIdx).';
A=[bLow -bLow.^2/6];
for i = 1:size(SigX,1)
    yi=-log(SigX(i,bIdx)).';
    DK=A\yi;
    DSigX(i,1)=DK(1);
    KSigX(i,1)=DK(2)/DK(1)^2;
    
    yi=-log(SigY(i,bIdx)).';
    DK=A\yi;
    DSigY(i,1)=DK(1);
    KSigY(i,1)=DK(2)/DK(1)^2;
end

% plot figure D and K of DKI

subplot(223);
hold on
hx=plot((TN*dt),DSigX,'bo'); axis square
set(hx,'markersize',6,'markerfacecolor',[0 0 1])
hy=plot((TN*dt),DSigY,'ro');
set(hy,'markersize',6)

% TD=dt*101:0.01:15;
% DTheory=1^2/4./TD;
% ht=plot(1./(TD),DTheory,'-k');
% set(ht,'linewidth',1);

% DMitra=d0*(1-4/6/sqrt(pi)*(2/rin)*sqrt(d0*TD));
htm=plot((TD),DMitra,'-.');
set(htm,'linewidth',1,'color',[0 0.6 0]);

hte=plot((TD),DTheoryExact,'--k');
set(hte,'linewidth',1);

legend([hx hy hte htm],{'$D_x({\rm DKI})$','$D_y({\rm DKI})$',...
    '$D_{\rm theory\,exact}$','$D_{\rm Mitra}$'},'interpreter','latex','fontsize',20,'location','southwest')

xlim([0 (max(TN)*dt)]); 
ylim([0 2.1]); box on;
xlabel('${\Delta}$ (ms)','interpreter','latex','fontsize',20);
ylabel('$D_\perp$ ($\mu$m$^2$/ms)','interpreter','latex','fontsize',20);
title(['$b=0.1-' num2str(bMax) '$ ms/$\mu$m$^2$'],'interpreter','latex','fontsize',20)
set(gca,'xscale','log')

subplot(224);
hold on
hx=plot(sqrt(TN*dt),KSigX,'bo'); axis square
set(hx,'markersize',6,'markerfacecolor',[0 0 1])
hy=plot(sqrt(TN*dt),KSigY,'ro');
set(hy,'markersize',6)
% ht=refline(0,-0.5);
% set(ht,'color',[0 0 0])
hte=plot(sqrt(TD),KTheoryExact,'--k');
set(hte,'linewidth',1);
% htm=plot(sqrt(TD),KMitra,'-.');
% set(htm,'linewidth',1,'color',[0 0.6 0])

legend([hx hy hte],{'$K_x({\rm DKI})$','$K_y({\rm DKI})$',...
    '$K_{\rm theory\,exact}$'},...
    'location','southeast','interpreter','latex','fontsize',20)

xlim([0 sqrt(max(TN)*dt)]); 
ylim([-2 3]); 
box on;
xlabel('$\sqrt{\Delta}$ (ms)','interpreter','latex','fontsize',20);
ylabel('$K_\perp$','interpreter','latex','fontsize',20);
title(['$b=0.1-' num2str(bMax) '$ ms/$\mu$m$^2$'],'interpreter','latex','fontsize',20)

%% plot figure Cumulant vs DKI

figure;
subplot(121);
hold on
hx=plot(DCumX,DSigX,'bo'); axis equal;
set(hx,'markerfacecolor',[0 0 1],'markersize',6);
hy=plot(DCumY,DSigY,'ro'); axis equal;
set(hy,'markersize',6);

xlabel('$D_{\rm Cumulant}$ ($\mu$m$^2$/ms)','interpreter','latex','fontsize',20)
ylabel('$D_{\rm DKI}$ ($\mu$m$^2$/ms)','interpreter','latex','fontsize',20)
hline=refline(1,0);
set(hline,'color',[0 0 0])
xlim([0 1.5]); ylim([0 1.5]);
box on;
legend([hx hy],{'$D_x$','$D_y$'},'location','southeast','interpreter','latex','fontsize',20)

subplot(122);
hold on
hx=plot(KCumX,KSigX,'bo'); axis equal;
set(hx,'markerfacecolor',[0 0 1],'markersize',6);
hy=plot(KCumY,KSigY,'ro'); axis equal;
set(hy,'markersize',6);

xlabel('$K_{\rm Cumulant}$','interpreter','latex','fontsize',20)
ylabel('$K_{\rm DKI}$','interpreter','latex','fontsize',20)
xlim([-1 3]); ylim([-1 3]);
hline=refline(1,0);
set(hline,'color',[0 0 0])
box on;
legend([hx hy],{'$K_x$','$K_y$'},'location','southeast','interpreter','latex','fontsize',20)

%%
Ntot=50:50:300;
% for i = 1:numel(Ntot)
N=700;%Ntot(i);
dJ0 = @(x) -besselj(1,x);
beta0k = @(k) fzero(dJ0, [(k-1)*pi, k*pi]); 
beta0=zeros(1,N); for k=1:N, beta0(k) = beta0k(k+1); end;
% beta0=[3.8317    7.0156   10.1735   13.3237   16.4706   19.6159   22.7601   25.9037   29.0468   32.1897];
% beta0 = [NaN beta0(1:N-1)];

% N=10;
dJ1 = @(x) besselj(0,x) - besselj(1,x)./(x + 1e-15);
beta1k = @(k) fzero(dJ1, [(k-1)*pi, k*pi]); 
beta1=zeros(1,N); for k=1:N, beta1(k) = beta1k(k); end;
% beta1 = [1.8412    5.3314    8.5363   11.7060   14.8636   18.0155   21.1644   24.3113   27.4571   30.6019];

% N=10;
dJ2 = @(x) besselj(1,x) - besselj(2,x)*2./(x + 1e-5);
beta2k = @(k) fzero(dJ2,[floor((k-1)*pi), k*pi]);
beta2=zeros(1,N); for k=1:N, beta2(k) = beta2k(k+1); end;
% beta2=[3.0542    6.7061    9.9695   13.1704   16.3475   19.5129   22.6716   25.8260   28.9777   32.1273];
% beta2 = [NaN beta2(1:N-1)];

% numer1=5/64/3+sum(1./beta0.^4+1./beta1.^2./(beta1.^2-1).*(4./beta1.^2-3/2)+1/2./beta2.^2./(beta2.^2-4));
% numer2=sum(-1./beta0.^2-1./(beta1.^2-1).*(4./beta1.^2-3/2)-1./(beta2.^2-4)/2);
% numer3=sum(1/2+1-3/4*beta1.^2./(beta1.^2-1)+1/4*beta2.^2./(beta2.^2-4))

ZeOr=5/64/3+sum(1./beta0.^4+1./beta1.^2./(beta1.^2-1).*(4./beta1.^2)+1./beta2.^2./(beta2.^2-4)/2)-3/16

frOr=sum(-1./beta0.^2-1/2./(beta2.^2-4))+1/4

% SeOr=
% end

%%

dis=res;
Nax=numel(r);

for i = 1:Nax-1
    for j = i+1:Nax
        temp=sqrt( (x(i)*res-x(j)*res)^2 + (y(i)*res-y(j)*res)^2 )-r(i)*res-r(j)*res;
        if temp<dis
            dis = temp;
            id=i;
            jd=j;
        end
%         display(['i=' num2str(i) ', j=' num2str(j)])
    end
end

%%













