function [x, y, rs] = initposition(rinit) 
%INITPOSITION    Initialize positions and radii of densly packed cylinders
%   [x,y,rs] = initposition(rinit) returns initial positions x, y and
%   rescaled radii rs of cylinders for Donev's C++ input file, based on
%   initial radii rinit (column vector), and saves them in the file
%   specified in readfilename.
%
%   ---------------------------------------------------------
%     in box.C, in the function
%     void box::ReadPositions(const char* filename)
%
%   infile.ignore(256, '\n');  // ignore the dim line
%   infile.ignore(256, '\n');  // ignore the #sphere 1 line
%   infile.ignore(256, '\n');  // ignore the #sphere line
%   infile.ignore(256, '\n');  // ignore the diameter line
%   infile.ignore(1000, '\n'); // ignore the 100 010 001 line
%   infile.ignore(256, '\n');  // ignore the T T T line
% 
%   for (int i=0; i<N; i++)
%     {
%       infile >> s[i].r;      // read in radius    
%       infile >> s[i].gr;     // read in growth rate
%       infile >> s[i].m;      // read in mass
%       for (int k=0; k<DIM; k++)  
%          infile >> s[i].x[k]; // read in position 
%     }
%    ..... 


readfilename='spheres_poly/read.dat';
dim=2;              % dimension
N=length(rinit);    % # cylinders

% rescale radii to make them small enough to not overlap, but not too small
rmax=max(rinit); 
dens0=0.01;         % dens0 ~ N(rmax/Rscale)^2
% divide all radii by Rscale for less initial overlap
Rscale=rmax*sqrt(N/dens0); 
rinit=sort(rinit(:),'descend');
rs=rinit/Rscale;

% assign random positions and check for no overlap
x=zeros(N,1); y=zeros(N,1); 
% distance^2 function
dist2 = @(x1,y1, x2,y2) min(abs(x1-x2),1-abs(x1-x2))^2 + min(abs(y1-y2),1-abs(y1-y2))^2;
for ncurr=1:N
   overlap=1;
   while (overlap==1)
      overlap=0;
      x(ncurr)=rand; y(ncurr)=rand;
      for nprev=1:ncurr-1
         if dist2(x(ncurr),y(ncurr),x(nprev),y(nprev))<=(rs(ncurr)+rs(nprev))^2, overlap=1; ncurr; break, end
      end
   end
end

% create read.dat file
fid=fopen(readfilename,'w');
fprintf(fid,'%d\n', dim);
fprintf(fid,'%d\n', N);
fprintf(fid,'%d\n', N);
fprintf(fid,'%e\n', 2*max(rs));
fprintf(fid,'10 01\n');
fprintf(fid,'T T\n');
for n=1:N
   fprintf(fid, '%e %e %f %f %f\n', rs(n), rs(n), 1.0, x(n), y(n));
end
fclose(fid);

end



