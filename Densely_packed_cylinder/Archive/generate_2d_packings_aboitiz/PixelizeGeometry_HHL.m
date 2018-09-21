function [A,B,res,Nmax]=PixelizeGeometry_HHL(res, xm, ym, r)
% function [A,B,res,Nmax]=PixelizeGeometry_HHL(res, xm, ym, r)
% Input:
%     res: about 1/( voxel# along one side )
%     xm,ym,r: positions and radii of axons
% Output:
%     A: axon index in voxels
%     B: axon # in voxels, either 1 or 2
%     res: as the input
%     Nmax: the smallest number, which is > total axon #, in the base of ten.
tic
% res = L/MS;
Xmax = 1; Ymax = 1;
MSX = ceil(Xmax/res)+1;
MSY = ceil(Ymax/res)+1;
A = zeros(MSX,MSY);
B = zeros(MSX,MSY);
surface = 0; %S2 = 0;
% h = waitbar(0,'pixelizing');
xmp=xm;
ymp=ym;
rp=r;
N=max(size(r));
Nmax=10^(ceil(log(N)/log(10)));

for i = 1:size(r,1)
%     waitbar(i/size(r,2));
    for ii = floor((xmp(i)-rp(i))/res+1):ceil((xmp(i)+rp(i))/res+1)
        if (ii>MSX)
           ti=ii-MSX;
        elseif (ii < 1)
            ti=ii+MSX;
        else
           ti=ii;
        end
        for jj = floor((ymp(i)-rp(i))/res+1):ceil((ymp(i)+rp(i))/res+1)
            if jj>MSY
                tj=jj-MSY;
            elseif (jj < 1)
                tj = jj+MSY;
            else
                tj=jj;
            end
%             if ((((ii*res-res/2-xmp(i))^2+(jj*res-res/2-ymp(i))^2) <= (rp(i))^2))
            if inside_circle(ii,jj,xmp(i),ymp(i),rp(i),res)
                if ~A(ti,tj)
                    A(ti,tj) = (i);
                    B(ti,tj) = (1);
                else
                    if ~floor(A(ti,tj)/Nmax) %1 other circle
                        A(ti,tj) = A(ti,tj)*Nmax +(i);
                        B(ti,tj) = 2;
                    else %if ~floor(A(ti+1,tj+1)/Nmax/Nmax) %2 other circles
                        display('more than 2 disks in one voxel');
                    end
                end                    
            end
        end
    end
    %S1 = S1 + (abs(round(xmp(i)+rp(i))-round(xmp(i)-rp(i)))+ abs(round(ymp(i)+rp(i))-round(ymp(i)-rp(i))))*2;
    surface = surface + 2*pi*r(i);
%     hold on; imagesc(A);
end
% close(h);
A=uint32(A);
B=uint16(B);
fprintf(' * Matrix filling done ! *\n');
fprintf(' ----------------------------\n');
toc
end

function inside = inside_circle(ii,jj,xmp,ymp,rp,res)
% function inside = inside_circle(ii,jj,xmp,ymp,rp,res)
% If the voxel (ii,jj) overlaps with the circle (xmp,ymp,rp), inside = 1,
% else, inside = 0;
% res: resolution
    v1 = (((ii*res-xmp)^2+(jj*res-ymp)^2) <= (rp)^2);
    v2 = (((ii*res-res-xmp)^2+(jj*res-ymp)^2) <= (rp)^2);
    v3 = (((ii*res-xmp)^2+(jj*res-res-ymp)^2) <= (rp)^2);
    v4 = (((ii*res-res-xmp)^2+(jj*res-res-ymp)^2) <= (rp)^2);
    
    inside = v1+v2+v3+v4 >0;
end

% figure; imagesc(A,[0 Nmax]); axis equal