function [A,B,Nmax] = createlookuptable(n,xc,yc,r)
%CREATELOOKUPTABLE    create lookup table for packed cylinders
%   [A,B,Nmax] = createlookuptable(n,xc,yc,r)
%   
%   Input:
%   n: size of the lookup table = n x n
%   xc, yc: center of cylinders, 0 <= xc,yc <= 1
%   r: radius of cylinders
%
%   Output:
%   A: axon labels/lookup table
%   B: # axon in each pixel
%   Nmax: the smallest integer larger than # axon, in the base of 10
%
% (c) Hong-Hsi Lee, 2016

MSX = n;
MSY = n;
A = zeros(MSX,MSY);
B = zeros(MSX,MSY);
N = length(r);
Nmax = 10^(ceil(log10(N)));

for i = 1:size(r,1)
    for ii = ceil((xc(i)-r(i))*n):ceil((xc(i)+r(i))*n)
        if (ii>MSX)
           ti = ii-MSX;
        elseif (ii < 1)
           ti = ii+MSX;
        else
           ti = ii;
        end
        for jj = ceil((yc(i)-r(i))*n):ceil((yc(i)+r(i))*n)
            if jj > MSY
                tj = jj-MSY;
            elseif (jj < 1)
                tj = jj+MSY;
            else
                tj = jj;
            end
            if inside_circle(ii,jj,xc(i),yc(i),r(i),n)
                if A(ti,tj) == 0        % 1 circle
                    A(ti,tj) = i;
                    B(ti,tj) = 1;
                else
                    if A(ti,tj) < Nmax  % 2 circles
                        A(ti,tj) = A(ti,tj)*Nmax + i;
                        B(ti,tj) = 2;
                    else                % > 2 circles
                        fprintf('More than 2 circles in one pixel.\n');
                    end
                end                    
            end
        end
    end
end
A = uint32(A);
B = uint16(B);
fprintf(' * Matrix filling done ! *\n');
fprintf(' ----------------------------\n');
end

function inside = inside_circle(ii,jj,xc,yc,r,n)
%INSIDE_CIRCLE    True for the pixel overlapping the circle
%   inside_circle(ii,jj,xc,yc,r,n) returns 1 if the pixel (ii,jj) overlaps
%   with the circle (xc,yc,r), otherwise, return 0. The size of the lookup
%   table is n x n.
%
% (c) Hong-Hsi Lee, 2016

    v1 = ((ii-xc*n)^2+(jj-yc*n)^2) <= (r*n)^2;
    v2 = ((ii-1-xc*n)^2+(jj-yc*n)^2) <= (r*n)^2;
    v3 = ((ii-xc*n)^2+(jj-1-yc*n)^2) <= (r*n)^2;
    v4 = ((ii-1-xc*n)^2+(jj-1-yc*n)^2) <= (r*n)^2;
    
    inside = ( (v1+v2+v3+v4) >0 );
end
