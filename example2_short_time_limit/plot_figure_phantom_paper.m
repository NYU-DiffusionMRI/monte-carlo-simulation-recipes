%% Plot Mitra limit of a randomly packed circles with Aboitiz's radius distributtion
close all

figure;
hold on
hx=plot(sqrt(TN*dt),(DCumX+DCumY)/2/d0,'-'); axis square
set(hx,'linewidth',3)

S = res*2*pi*sum(r);
V = res^2*(1-pi*sum(r.^2));
SoV= S/V;
DMitra = 1-SoV*2/3*sqrt(d0*TN*dt)/sqrt(pi);
KMitra=3/5*SoV*sqrt(d0*TN*dt)/sqrt(pi);

hdm = plot(sqrt(TN*dt),DMitra','--r'); set(hdm,'linewidth',2);
legend([hx,hdm],{'Simulation','Mitra limit'},'interpreter','latex','fontsize',30,'location','southwest')

set(gca,'xtick',0:0.05:1,'ytick',0:0.5:1,'fontsize',20);
xlim([0 0.3]); 
ylim([0 1.25]); box on; grid on
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',30);
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',30);

figure;
hold on
hx=plot(sqrt(TN*dt),(KCumX+KCumY)/2,'-'); axis square
set(hx,'linewidth',3)
hkm=plot(sqrt(TN*dt),KMitra,'--r');
set(hkm,'linewidth',2);
hl = refline(0,0); set(hl,'color','k')
legend([hx,hkm],{'Simulation','Mitra limit'},'interpreter','latex','fontsize',30,'location','southeast')

set(gca,'xtick',0:0.05:1,'ytick',-2:0.5:1,'fontsize',20);
xlim([0 0.3]); 
ylim([-2 1]); box on; grid on
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',30);
ylabel('$K(t)$','interpreter','latex','fontsize',30);

%% Plot geometry

figure;
for i = -1:1
    for j = -1:1
        viscircles([x(:)+i,y(:)+j],r(:),'color','k');
    end
end
xlim([0 1]);
ylim([0 1]);
set(gca,'xtick',[],'ytick',[],'linewidth',3);
box on
axis square

figure;
viscircles([0.5,0.5],0.28,'color','k');
xlim([0 1]);
ylim([0 1]);
set(gca,'xtick',[],'ytick',[],'linewidth',3);
box on
axis square
%%
addpath('/Users/magda/Desktop/diffusion_simulation/cpp_practice/generate_packing')
[A,B,Np,Nmax] = PixelizeGeometry_HHL_v2(4e2, x, y, r);

cmap = colormap('parula');
Nax  = 492;
Ibg = A==0;
Iol = A>Nax;
A2 = ceil(single(A)/492*64);
A2(Ibg) = 1;
A2(Iol) = 1;

[nx,ny] = size(A2);
imgr = reshape(A2,[nx*ny,1]); imgg = imgr; imgb = imgr;
imgr = cmap(uint16(imgr),1); imgr = reshape(imgr,[nx,ny]);
imgg = cmap(uint16(imgg),2); imgg = reshape(imgg,[nx,ny]);
imgb = cmap(uint16(imgb),3); imgb = reshape(imgb,[nx,ny]);

imgr(Ibg) = 0; imgg(Ibg) = 0; imgb(Ibg) = 0;
imgr(Iol) = 1; imgg(Iol) = 1; imgb(Iol) = 1;
imgc = cat(3,imgr,imgg,imgb);

figure
image(rot90(imgc));
box on; axis square
set(gca,'xtick',[],'ytick',[],'linewidth',3)
colormap parula
h = colorbar; set(h,'fontsize',12); caxis([0 492])





