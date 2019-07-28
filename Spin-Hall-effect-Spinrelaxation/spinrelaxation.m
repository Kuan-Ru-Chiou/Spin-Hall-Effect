%spin-hall-effect-rashba prl paper parameter

clear all;close all;

%constant
 hbar = 1.06e-34;    %plank constant
 q = 1.6e-19;       %fundametal electron charge
 m = 0.05*9.1e-31;  %electron effective mass
 qh = q/hbar;
zplus = 1i*1e-15;    %green function use small imaginary part
kT = 1.d-10;         %thermal energy
% DD = 0;       %momentum or phase relaxation strength
Ds =0.1;    %spin-relaxation
%pauli matrix
sx = [0 1; 1 0];
sy = [0 -1i;1i 0];
sz = [1 0;0 -1];

%input 
a = 3e-9; %lattice constant
t0 = (hbar^2)/(2*m*(a^2)*q); % hopping energy
eta = a*0.1*t0; %rashba strength   
NW = 11; % dervice width length
%n =14;
Np =11;  %device transport length
NWp = NW*Np; %total grid points
Ef = -3.8*t0; %fermi energy

%%-------------------------

%%initialize left and right contact self energy
L = zeros(Np);
R = L;
L(1,1) = 1;
R(Np,Np) = 1;


%%% construct  hamitonian
al = 0*t0*eye(2); by = -t0*eye(2) - (1i*eta/2/a)*sx;    % band bottom 
byy = -t0*eye(2);
bx = -t0*eye(2)+(1i*eta/2/a)*sy;                        
bxx = -t0*eye(2);

alphaa = kron(eye(NW),al) +kron(diag(ones(1,NW-1),1),byy)+kron(diag(ones(1,NW-1),-1),byy'); %%for contact no rashba
alpha = kron(eye(NW),al) +kron(diag(ones(1,NW-1),1),by)+kron(diag(ones(1,NW-1),-1),by');  %%% for  contact contain rashba
beta = kron(diag(ones(1,NW)),bx);    %%% for  contact contain rashba
betaa = kron(diag(ones(1,NW)),bxx);  %%for contact no rashba

H = kron(eye(Np),alpha);

if Np>1 
    H = H+kron(diag(ones(1,Np-1),1),beta)+kron(diag(ones(1,Np-1),-1),beta');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%TB band structure
ka = linspace(-pi,pi,501);

for k = 1:length(ka)
    [V,D] = eig(alpha+beta*exp(1i*ka(k))+beta'*exp(-1i*ka(k)));
    eigE(:,k) = sort(diag(D));
end
figure(7);
plot(ka/pi,eigE);




ii = 1;
A0 = zeros(NWp,1);
Ax = zeros(2*NWp,2*NWp);
Ay = zeros(2*NWp,2*NWp);
Az = zeros(2*NWp,2*NWp);


%for spin relaxation use
Sxxx = kron(eye(NWp),sx);
Syyy = kron(eye(NWp),sy);
Szzz = kron(eye(NWp),sz);



%bias
V = 1.d-3*t0;    %if want to test oppisite bias change to mins

mu1 = Ef + V/2; mu2 = Ef -V/2;

%NEGF calculation for every energy
for EE = linspace(mu2,mu1,1)
    
% sigB =zeros(2*NWp); sigBin = sigB; %initialize guess value for dephasing
    
 sigBs = zeros(2*NWp); sigBsin = sigBs; 
    
%  f1=fdis(EE,mu1,kT);
%  f2=fdis(EE,mu2,kT);
f1 =1;  % if want to test oppisite bias change to 0 and V change to -V
f2 =0;   % if want to test oppisite bias change to 1 and V change to -V


 %%%%%%Sancho-ruby iteration solve surface green function 
   
     g1alpha = (EE+zplus)*eye(2*NW) - alphaa - diag(ones(1,2*NW)*(V/2));
     g1beta = -(betaa)';
     change = 1;
     t = inv(g1alpha)*(betaa);  % initial value
     td = inv(g1alpha)*(betaa)';% initial value 
     T = t;
     Mt = td;
     %criterion is all the renormalization hopping is smaller than change
    while change > 1e-9
        tt = inv(eye(2*NW) - t*td-td*t)*t*t;    
        ttd = inv(eye(2*NW) - t*td - td*t)*td*td; 
        S = T + Mt*tt;
        change = sum(sum(abs(S-T)))/sum(sum(abs(S+T)));
        T = S;
        t = tt;
        td = ttd;
        Mt = Mt*ttd;
    end
    g1 = inv(g1alpha - (betaa')*T);
    g1 = (betaa)*g1*(betaa)';
    sig1 = kron(diag([1 zeros(1,Np-1)]),g1);
    gam1 = 1i*(sig1 - sig1');
    sig1in = gam1*f1;
    
    g2alpha = (EE+zplus)*eye(2*NW) - alphaa - diag(ones(1,2*NW)*(-V/2));
     g2beta = -(betaa);
     change = 1;
     t = inv(g2alpha)*(betaa');
     td = inv(g2alpha)*(betaa);
     T = t;
     Mt = td;
     
    while change > 1e-9
        tt = inv(eye(2*NW) - t*td-td*t)*t*t;
        ttd = inv(eye(2*NW) - t*td - td*t)*td*td;
        S = T + Mt*tt;
        change = sum(sum(abs(S-T)))/sum(sum(abs(S+T)));
        T = S;
        t = tt;
        td = ttd;
        Mt = Mt*ttd;
    end
    g2 = inv(g2alpha - (betaa)*T);
    g2 = (betaa)*g2*(betaa)';
    sig2 = kron(diag([zeros(1,Np-1) 1]),g2);
    gam2 = 1i*(sig2 - sig2');
    sig2in = gam2*f2;
 
       %green function self consistent for dephasing (momentum or phase)
        
%      change = 100;
%      while change > 1e-6
% 
%          G = inv((EE*eye(2*NWp))-H-sig1-sig2-sigB);
%          sigBnew = DD*eye(2*NWp).*G;
%          change = sum(sum(abs(sigB-sigBnew)))/sum(sum(abs(sigB+sigBnew)));
%          sigB = sigB + 0.25*(sigBnew-sigB);
%      end
%      
%      
%        change = 100;
%      while change > 1e-6
%          gn = G*(sig1in+sig2in+sigBin)*G';
%          sigBinnew = DD*eye(2*NWp).*gn;
%         change = sum(sum(abs(sigBin-sigBinnew)))/sum(sum(abs(sigBin+sigBinnew)));
%          sigBin = sigBin +0.25*(sigBinnew - sigBin) ;
%      end
  


    %green function self consistent for dephasing (spin relaxation) 
     GG = zeros(2*NWp); ggn =GG;
     change = 100;
     while change > 1e-12

         G = inv((EE+zplus)*eye(2*NWp)-H-sig1-sig2-sigBs);
         for i =1:NWp
             GG(2*i-1:2*i,2*i-1:2*i) = G(2*i-1:2*i,2*i-1:2*i);
         end
         sigBsnew = Ds*(Sxxx*GG*Sxxx +Syyy*GG*Syyy+Szzz*GG*Szzz);
         change = sum(sum(abs(sigBs-sigBsnew)))/sum(sum(abs(sigBs+sigBsnew)));
         sigBs = sigBs + 0.5*(sigBsnew-sigBs);
     end
     
     
       change = 100;
     while change > 1e-12
         gn = G*(sig1in+sig2in+sigBsin)*G';
         for i =1:NWp
             ggn(2*i-1:2*i,2*i-1:2*i) = gn(2*i-1:2*i,2*i-1:2*i);
         end
         sigBsinnew = Ds*(Sxxx*ggn*Sxxx +Syyy*ggn*Syyy+Szzz*ggn*Szzz);
        change = sum(sum(abs(sigBsin-sigBsinnew)))/sum(sum(abs(sigBsin+sigBsinnew)))
         sigBsin = sigBsin +0.5*(sigBsinnew - sigBsin) ;
     end
     
     
     
     
     


 
 %%%%using keldish Gless formula to check

Gless =  1i*gn;
%%%%%%%

 
 A = 1i*(G-G');
%   Gn = G*(gam1in+gam2in)*G';
 Tcosh(ii) = real(trace(gam1*G*gam2*G'));


 
ii = ii+1;

%   E(ii) = EE/t0;
end
%%%%%%%%%%%%%%%%%%%%%%calculation spin density
N0=zeros(NWp,1);
Nx=zeros(NWp,1);
Ny=zeros(NWp,1);
Nz=zeros(NWp,1);
A00 = zeros(NWp,1);

for i=1:length(N0)
    N0(i,1) = (-1i)*trace(Gless((2*i-1):2*i,(2*i-1):2*i)*eye(2,2));
end
for i=1:length(Nx)
    Nx(i,1) = (-1i)*trace(Gless((2*i-1):2*i,(2*i-1):2*i)*sx);
end
for i=1:length(Ny)
    Ny(i,1) =(-1i)*trace(Gless((2*i-1):2*i,(2*i-1):2*i)*sy);
end
for i=1:length(Nz)
    Nz(i,1) = (-1i)*trace(Gless((2*i-1):2*i,(2*i-1):2*i)*sz);
end
for i =1:length(A00)
 A00(i,1) = trace(A((2*i-1):2*i,(2*i-1):2*i));
end

 
 fx = Nx./A00;
 fy = Ny./A00;
 fz = Nz./A00;
 f0 = N0./A00;

 
  Y = [0:1:NW-1]*a;
%check mid-point chanel quasi-fermi-level transport
  f00 =reshape(f0,NW,NW);
  f00 = f00(round(NW/2),:);
  
  fxx =reshape(fx,NW,NW);
  fxx = fxx(round(NW/2),:);
  
  fyy =reshape(fy,NW,NW);
  fyy = fyy(round(NW/2),:);
  
  fzz =reshape(fz,NW,NW);
  fzz = fzz(round(NW/2),:);
  
  figure(14)
  plot(Y,f00);
  hold on
  plot(Y,fxx);
  hold on
  plot(Y,fyy);
  hold on
  plot(Y,fzz);
  hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

 
 
 
 
 
 
 
 



%%%%%%%%%%%%%%%%%%%%%%



%%%%%  contour  plot for spin density
[XX,YY] = meshgrid([0:1:Np-1],[0:1:NW-1]);
XX=a*XX;
YY=a*YY;
Y = [0:1:NW-1]*a;





F0 = real(N0);
Fx = real(Nx);
Fy = real(Ny);
Fz = real(Nz);


%check charge conservation in transport cross section

for i=1:NWp-NW
   J0(i) = (-1/2)*trace(H(2*(i+NW)-1:2*(i+NW),2*i-1:2*i)*Gless(2*i-1:2*i,2*(i+NW)-1:2*(i+NW)) - H(2*i-1:2*i,2*(i+NW)-1:2*(i+NW))*Gless(2*(i+NW)-1:2*(i+NW),2*i-1:2*i));
end
J0 = real(J0);
J0 = reshape(J0,NW,Np-1);
[HH,RR] = meshgrid([1:1:Np-1],[1:1:NW]);
quiver(HH,RR,J0,zeros(NW,Np-1));
title('charge current');
fprintf(2,'check charge current = %g\n',sum(J0,1));  %check current charge current conservation


%%%%%%%%%%%check terminal charge and spin current
JCC1= trace(sig1in*A-gam1*gn);
JCC2= trace(sig2in*A-gam2*gn);

JJZ1 =trace((sig1in*A-gam1*gn)*Szzz);
JJZ2 =trace((sig2in*A-gam2*gn)*Szzz);

JJY1 =trace((sig1in*A-gam1*gn)*Syyy);
JJY2 =trace((sig2in*A-gam2*gn)*Syyy);

JJX1 = trace((sig1in*A-gam1*gn)*Sxxx);
JJX2 = trace((sig2in*A-gam2*gn)*Sxxx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check spin-current in transport cross section must be conservation
for i=1:NWp-NW
   Jzz0(i) = (-1/2)*trace((H(2*(i+NW)-1:2*(i+NW),2*i-1:2*i)*sz + sz*H(2*(i+NW)-1:2*(i+NW),2*i-1:2*i))*Gless(2*i-1:2*i,2*(i+NW)-1:2*(i+NW)) - (H(2*i-1:2*i,2*(i+NW)-1:2*(i+NW))*sz + sz*H(2*i-1:2*i,2*(i+NW)-1:2*(i+NW)))*Gless(2*(i+NW)-1:2*(i+NW),2*i-1:2*i));
end

for i=1:NWp-NW
   Jyy0(i) = (-1/2)*trace((H(2*(i+NW)-1:2*(i+NW),2*i-1:2*i)*sy + sy*H(2*(i+NW)-1:2*(i+NW),2*i-1:2*i))*Gless(2*i-1:2*i,2*(i+NW)-1:2*(i+NW)) - (H(2*i-1:2*i,2*(i+NW)-1:2*(i+NW))*sy + sy*H(2*i-1:2*i,2*(i+NW)-1:2*(i+NW)))*Gless(2*(i+NW)-1:2*(i+NW),2*i-1:2*i));
end

for i=1:NWp-NW
   Jxx0(i) = (-1/2)*trace((H(2*(i+NW)-1:2*(i+NW),2*i-1:2*i)*sx + sx*H(2*(i+NW)-1:2*(i+NW),2*i-1:2*i))*Gless(2*i-1:2*i,2*(i+NW)-1:2*(i+NW)) - (H(2*i-1:2*i,2*(i+NW)-1:2*(i+NW))*sx + sx*H(2*i-1:2*i,2*(i+NW)-1:2*(i+NW)))*Gless(2*(i+NW)-1:2*(i+NW),2*i-1:2*i));
end


Jzz0 = real(Jzz0);
Jzz0 = reshape(Jzz0,NW,Np-1);

Jxx0 = real(Jxx0);
Jxx0 = reshape(Jxx0,NW,Np-1);

Jyy0 = real(Jyy0);
Jyy0 = reshape(Jyy0,NW,Np-1);

[HH,RR] = meshgrid([1:1:Np-1],[1:1:NW]);
quiver(HH,RR,Jzz0,zeros(NW,Np-1));
title('spin current sz');


[HH,RR] = meshgrid([1:1:Np-1],[1:1:NW]);
quiver(HH,RR,Jxx0,zeros(NW,Np-1));
title('spin current sx');


[HH,RR] = meshgrid([1:1:Np-1],[1:1:NW]);
quiver(HH,RR,Jyy0,zeros(NW,Np-1));
title('spin current sy');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%spin current in transverse direction
Jyyy0 = zeros(1,NW*Np); Jzyy0=Jyyy0; Jxyy0=Jyyy0 ;
for i = 1:NW*Np-1
    
    Jzyy0(i) = (-1/2)*trace((H(2*i+1:2*i+2,2*i-1:2*i)*sz + sz*H(2*i+1:2*i+2,2*i-1:2*i))*Gless(2*i-1:2*i,2*i+1:2*i+2) - (H(2*i-1:2*i,2*i+1:2*i+2)*sz + sz*H(2*i-1:2*i,2*i+1:2*i+2))*Gless(2*i+1:2*i+2,2*i-1:2*i)); 
    Jyyy0(i) = (-1/2)*trace((H(2*i+1:2*i+2,2*i-1:2*i)*sy + sy*H(2*i+1:2*i+2,2*i-1:2*i))*Gless(2*i-1:2*i,2*i+1:2*i+2) - (H(2*i-1:2*i,2*i+1:2*i+2)*sy + sy*H(2*i-1:2*i,2*i+1:2*i+2))*Gless(2*i+1:2*i+2,2*i-1:2*i)); 
    Jxyy0(i) = (-1/2)*trace((H(2*i+1:2*i+2,2*i-1:2*i)*sx + sx*H(2*i+1:2*i+2,2*i-1:2*i))*Gless(2*i-1:2*i,2*i+1:2*i+2) - (H(2*i-1:2*i,2*i+1:2*i+2)*sx + sx*H(2*i-1:2*i,2*i+1:2*i+2))*Gless(2*i+1:2*i+2,2*i-1:2*i)); 
end

Jyyy0(NW*Np) = 0; Jzyy0(NW*Np) = 0; Jxyy0(NW*Np) = 0;

Jzyy0 = real(Jzyy0);
Jzyy0 = reshape(Jzyy0,NW,Np);

Jyyy0 = real(Jyyy0);
Jyyy0 = reshape(Jyyy0,NW,Np);

Jxyy0 = real(Jxyy0);
Jxyy0 = reshape(Jxyy0,NW,Np);

%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%plot spin current quiver

Jzyy0(:,Np)=[];
Jyyy0(:,Np)=[];
Jxyy0(:,Np)=[];


figure(11);
quiver(HH,RR,Jzz0,Jzyy0);
title('spin current sz');
figure(12);
quiver(HH,RR,Jyy0,Jyyy0);
title('spin current sy');
figure(13);
quiver(HH,RR,Jxx0,Jxyy0);
title('spin current sx');



% I1 = real(trace(gam1*f1*A+gam1*1i*Gless));  %check terminal 1 current 
% I2 = -1*real(trace(gam2*f2*A+gam2*1i*Gless)); %check terminal 2 current  for contact 2 be careful minus sign
% fprintf(2,'check terminal 1 charge current = %g\n',I1); 
% fprintf(2,'check terminal 2 charge current = %g\n',I2); 





%%%%%%%%%%%plot spin density contour in 2D----------------------
% tempx=(1:Np)';
% x=kron(tempx,ones(NW,1));
% y=repmat((1:NW)',Np,1);
% z=Fz(:);
% 
% [XX YY ZZ]=gridplot(x,y,z,'x','y','spin z','contourf');
%%%%%%%%%%%%

Szz = reshape(Fz,NW,Np);
figure(1);
contourf(XX,YY,Szz,50,'linestyle','none'); shading flat; colormap(hot(50));
colorbar;
title('Szz spin accumulation');
xlabel('lontudinal (m)')
ylabel('transverse (m)')

Sxx = reshape(Fx,NW,Np);
figure(2);
contourf(XX,YY,Sxx,50,'linestyle','none'); shading flat;colormap(hot(50));
colorbar;
title('Sxx spin accumulation');
xlabel('lontudinal (m)')
ylabel('transverse (m)')

Syy = reshape(Fy,NW,Np);
figure(3);
contourf(XX,YY,Syy,50,'linestyle','none'); shading flat; colormap(hot(50));
colorbar;
title('Syy spin accumulation');
xlabel('lontudinal (m)')
ylabel('transverse (m)')


%%%%%%%%%plot spin density in trasverse direction-----------------
% figure(4);plot(Y,Fz([1:NW]),'o-r','linewidth',2);
% xlabel('width(m)');
% ylabel('Sz accumulation(1/eV)');

figure(4);plot(Y,Fz([(round(NW/2)-1)*NW+1:round(NW/2)*NW]),'o-r','linewidth',2);
xlabel('transverse (m)');
ylabel('Sz accumulation (1/eV)');



figure(5);
plot(Y,Fx([(round(NW/2)-1)*NW+1:round(NW/2)*NW]),'o-r','linewidth',2);
xlabel('transverse (m)');
ylabel('Sx accumulation (1/eV)');

figure(6);
plot(Y,Fy([(round(NW/2)-1)*NW+1:round(NW/2)*NW]),'o-r','linewidth',2);
xlabel('transverse (m)');
ylabel('Sy accumulation (1/eV)');



figure(7);
 plot(Y,fx([(round(NW/2)-1)*NW+1:round(NW/2)*NW]),'o-r','linewidth',2);
xlabel('transverse (m)');
ylabel('fsx Occupation');
figure(8);
 plot(Y,fy([(round(NW/2)-1)*NW+1:round(NW/2)*NW]),'o-r','linewidth',2);
xlabel('transverse (m)');
ylabel('fsy Occupation');

figure(9);
 plot(Y,fz([(round(NW/2)-1)*NW+1:round(NW/2)*NW]),'o-r','linewidth',2);
xlabel('transverse (m)');
ylabel('fsz Occupation');

figure(10);
plot(Y,f0([(round(NW/2)-1)*NW+1:round(NW/2)*NW]),'o-r','linewidth',2);
xlabel('transverse (m)');
ylabel('fs0 Occupation');





