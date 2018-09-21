close all
n = 10000;

th = 2*pi*rand([n 1]);
psi = 2*pi*rand([n 1]);

x = sin(th).*cos(psi);
y = sin(th).*sin(psi);
z = cos(th);

figure; 

subplot(241)
scatter3(x,y,z,0*x+30,'marker','.'); pbaspect([1 1 1])
title('random \theta, \psi')

subplot(242)
scatter(x,y,0*x+30,'marker','.'); pbaspect([1 1 1])
title('xy projection')

subplot(243)
scatter(y,z,0*y+30,'marker','.'); pbaspect([1 1 1])
title('yz projection')

rxy = sqrt(x.^2+y.^2);
ryz = sqrt(y.^2+z.^2);
binrange = [0.05:0.01:1-0.05].';
Nxy = histc(rxy,binrange);
Nyz = histc(ryz,binrange);

subplot(244)
hold on
hxy = plot(binrange,Nxy./binrange,'-b'); hyz = plot(binrange,Nyz./binrange,'-r');
legend([hxy hyz],{'xy','yz'})
pbaspect([1 1 1])
xlabel('velocity'); ylabel('count')


xi = rand([n 3])-0.5;
ri = sqrt(sum(xi.^2,2));
ri = repmat(ri,[1 3]);
xi = xi./ri;

x = xi(:,1); y = xi(:,2); z = xi(:,3);

subplot(245)
scatter3(x,y,z,0*x+30,'marker','.'); pbaspect([1 1 1])
title('random vx,vy,vz')

subplot(246)
scatter(x,y,0*x+30,'marker','.'); pbaspect([1 1 1])
title('xy projection')

subplot(247)
scatter(y,z,0*y+30,'marker','.'); pbaspect([1 1 1])
title('yz projection')

rxy = sqrt(x.^2+y.^2);
ryz = sqrt(y.^2+z.^2);
binrange = [0.05:0.01:1-0.05].';
Nxy = histc(rxy,binrange);
Nyz = histc(ryz,binrange);

subplot(248)
hold on
hxy = plot(binrange,Nxy./binrange,'-b'); hyz = plot(binrange,Nyz./binrange,'-r');
legend([hxy hyz],{'xy','yz'},'location','northwest')
pbaspect([1 1 1])
xlabel('velocity'); ylabel('count')

%%
xi = 2*rand([n 2])-1;
xi = xi./sqrt(2);

x = xi(:,1); y = xi(:,2);

figure;
subplot(231)
scatter3(x,y,z,0*x+30,'marker','.'); pbaspect([1 1 1])

subplot(232)
scatter(x,y,0*x+30,'marker','.'); pbaspect([1 1 1])

r = sqrt(x.^2+y.^2);
binrange = [0.05:0.01:1-0.05].';
N= histc(r,binrange);

subplot(233)
plot(binrange,N./binrange)
