function [A,B,Np,Nmax]=PixelizeGeometry_HHL_v2(Np, xm, ym, r)
% function [A,B,Np,Nmax]=PixelizeGeometry_HHL_v2(Np, xm, ym, r)
% Input:
%     Np: pixel # along each side
%     xm,ym,r: positions and radii of axons
% Output:
%     A: axon index in each pixel
%     B: axon # in each pixel, either 1 or 2
%     Np: as the input
%     Nmax: the smallest number, which is > total axon #, in the base 10
tic

MSX = Np;
MSY = Np;
A = zeros(MSX,MSY);
B = zeros(MSX,MSY);
xmp=xm;
ymp=ym;
rp=r;
N=max(size(r));
Nmax=10^(ceil(log(N)/log(10)));

for i = 1:size(r,1)
    for ii = ceil((xmp(i)-rp(i))*Np):ceil((xmp(i)+rp(i))*Np)
        if (ii>MSX)
           ti=ii-MSX;
        elseif (ii < 1)
           ti=ii+MSX;
        else
           ti=ii;
        end
        for jj = ceil((ymp(i)-rp(i))*Np):ceil((ymp(i)+rp(i))*Np)
            if jj>MSY
                tj = jj-MSY;
            elseif (jj < 1)
                tj = jj+MSY;
            else
                tj=jj;
            end
            if inside_circle(ii,jj,xmp(i),ymp(i),rp(i),Np)
                if A(ti,tj) == 0
                    A(ti,tj) = i;
                    B(ti,tj) = 1;
                else
                    if A(ti,tj) < Nmax %1 other circle
                        A(ti,tj) = A(ti,tj)*Nmax + i;
                        B(ti,tj) = 2;
                    else %2 other circles
                        display('more than 2 circles in one pixel');
                    end
                end                    
            end
        end
    end
end
A=uint32(A);
B=uint16(B);
fprintf(' * Matrix filling done ! *\n');
fprintf(' ----------------------------\n');
toc
end

function inside = inside_circle(ii,jj,xmp,ymp,rp,Np)
% function inside = inside_circle(ii,jj,xmp,ymp,rp,Np)
% If the voxel (ii,jj) overlaps with the circle (xmp,ymp,rp), inside = 1,
% else, inside = 0;
% Np: pixel # along each side
    v1 = (((ii-xmp*Np)^2+(jj-ymp*Np)^2) <= (rp*Np)^2);
    v2 = (((ii-1-xmp*Np)^2+(jj-ymp*Np)^2) <= (rp*Np)^2);
    v3 = (((ii-xmp*Np)^2+(jj-1-ymp*Np)^2) <= (rp*Np)^2);
    v4 = (((ii-1-xmp*Np)^2+(jj-1-ymp*Np)^2) <= (rp*Np)^2);
    
    inside = ( (v1+v2+v3+v4) >0 );
end

% figure; imagesc(A,[0 Nmax]); axis equal