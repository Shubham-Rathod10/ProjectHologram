%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%%%%%  SIMULATION : RANDOM PARTICLE OBJECT AT A SINGLE PLANES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T0, hologramO] = funcParticle(S,nP)
%S = input('Insert the number of Object = '); % no. of object
%nP = input('Insert the number of particles in a single plane = '); % no. of particles in objects
for ii = 1:S
N = 1000;  % number of pixels
lambda = 633*10^(-9); % wavelength in meter
area = 0.005; % area side length in meter
zp = 0.08; % z in meter
r =  0.015*ones(nP,1); % radius of point objects

T0 =0.8*ones(N,N);
xc = -1+2*1*rand(nP,1); 
yc = -1+2*1*rand(nP,1);
T0 = funcParticleSim(T0,xc,yc,r);
%figure; imagesc(abs(T0)); colormap gray; %title('Object at Plane A'); 
prop0 = PropagatorN(N,lambda,area,zp);

imagesc(real((prop0)));colormap gray
hol1 = IFT2Dc(FT2Dc(T0).*prop0);
U0 = abs(hol1).^2;
hologramO = U0;
%figure; imagesc(hologramO);colormap gray; title('Hologram');

%%%%%%%%%%% Save Output %%%%%%%%%%%%%%%%%%%%%%%%%%
NameO = strcat('Object',int2str(ii),'.png');
NameH = strcat('Hologram',int2str(ii),'.png');
imwrite (abs(T0),NameO); colormap gray;
imwrite (hologramO,NameH); colormap gray;
end
end
%%
%%%%%%%%%%%  SUPPORTIVE %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% PropagatorN %%%%%%%%%%%%%%%%%%%%%%%%
function[ p ] = PropagatorN(N,lambda,area,z)
p = zeros(N,N);
ii= 1:N;
alpha = lambda*(ii-N/2-1)/area;
jj= 1:N;
beta = lambda*(jj-N/2-1)/area;
[alpha,beta] = meshgrid(alpha,beta);
p((alpha^2 + beta^2)<=1) =  exp(-2*pi*1i*z*sqrt(1 - alpha.^2 - beta.^2)/lambda);
% for ii = 1:N
%     for jj = 1:N
%         alpha = lambda*( ii - N/2-1)/area;
%         beta = lambda*( jj - N/2-1)/area;
%         if((alpha^2 + beta^2)<=1)
%      % p( ii , jj ) = exp(-2*pi*1i*z*sqrt(1 - alpha^2 - beta^2)/lambda);
%         end  % if
%     end
% end
end
%%%%%%%%%%%  2D Fourier Transform %%%%%%%%%%%%%%%%%%%%
function [out] = FT2Dc(in)
[M,N] = size(in);
f1 = zeros(M,N);
for ii = 1:M
    for jj = 1:N
         f1(ii,jj) = exp(1i*pi*(ii + jj));
         %in(ii,jj);
         %f1;
    end
end
% size(f1)
% size(in)
f2=f1.*in;
FT = fft2(f2);
out = f1.*FT;
end
%%%%%%%%%%%%%%%%%% 2D Inverse Fourier Transform %%%%%%%%%%%%%
function [out] = IFT2Dc(in)
[M,N] = size(in);
f1 = zeros(M,N);
for ii = 1:M
    for jj = 1:N
        f1(ii, jj) = exp(-1i*pi*(ii + jj));
    end
end
FT = ifft2(f1.*in);
out = f1.*FT;
end
%%%%%%%%%%%%% RANDOM PARTICLE GENERATION %%%%%%%%%%%%%%%%%%%%%%%%
function object = funcParticleSim(object,xc,yc,r)
[M,N]=size(object);
x = linspace(-1,1,N);
y = linspace(-1,1,M);
[X,Y] = meshgrid(x,y);
K = numel(xc);
if numel(xc) ~= numel(yc)
   msg = 'x and y coordinate vector should be of same length';
   error(msg)
end
if numel(xc) ~= numel(r)
   msg = 'x,y coordinate vectors and r should be of same length';
   error(msg)
end
for k = 1:K
    object((X-xc(k)).^2+(Y-yc(k)).^2 < r(k)^2) = 0;
  %  phase((X-xc(k)).^2+(Y-yc(k)).^2 > r(k)^2) = 0;
   %r object = exp(1i*phase);
end

end
%%
