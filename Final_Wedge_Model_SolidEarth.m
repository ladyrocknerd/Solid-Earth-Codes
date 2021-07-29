%% Wedge Kinematics

clear all
close all
%Given a wedge with a frontal flux in, underplating, erosion and rearward flux out. Assume steady
%state, so that the fluxes balance.

cd 'C:\Users\erlanger\Desktop\Matlab\Wedge_Model'
cd 'C:\Users\erlanger\Desktop\Matlab'

%%  Give following parameters:
h0=10.;  % incoming section for frontal accretion In %km
h1=10.;  % incoming section for prowedge basal accretion
h2=0.;  % incoming section for retrowedge basal accretion
ep=0.75;  % pro-wedge erosion rate%km/Ma 
er=0.25;  % retro-wedge erosion rate 
LP=60.; % length of Pro-wedge % km
LR=40.;  % Length of retro-Wedge
vp=8;  %Slab retreat rate %km/Ma  (Normally signifies convergence rate but Apennines has negligible convergence)
hm=56.; % Max Crustal Thickness
elevm=2;  % Maximum Elevation (km)

%% Wedge geometry

%prowedge slopes
alphap=elevm/LP;
betap=(hm-elevm-h0)/LP;
thetap=alphap+betap;

%determine thickness of outgoing crust hR and retro- slope
hR=h0+h1+h2-ep*LP/vp-er*LR/vp;
retrocrust=hR;

% Retro-wedge
alphar=(elevm)/LR;
betar=(hm-elevm-hR)/LR;
thetar=alphar+betar;

%alphap=tan(2*pi/180);
%betap=tan(12*pi/180);
%thetar=tan(8*pi/180);

x=linspace(0,LP,1001);
v=zeros(1,1001);
u=zeros(1001,101);
xr=linspace(0,LR,1001);

% Underplating velocity on prowedge
UP=vp*h1/LP;

%Underplating velocity on retrowedge
UR=vp*h2/LR;

%determine maximum thickness

hm=h0+LP*(alphap+betap);
MaxCrustThick=hm
elevm=LP*(alphap)

z=linspace(0,hm,101);

%% prowedge field

%find horizontal velocity across pro-wedge
v=(h0*vp+(UP-ep)*x)./(h0+x*thetap);

%find vertical velocity on the prowedge surface (defined as the erosion velocity plus the motion to carry the point
%along the surface)
us=ep+((h0*vp+x*(UP-ep))*alphap)./((h0+x*thetap));

% zs=h0+x*alphap; why are there two calculations of this?

% elevation of prowedge surface 
zs=x*alphap; 

% elevation of the base of the prowedge
zb=-h0-x*betap;

%% retrowedge field

% velocity at the divide
vm=(h0*vp+(UP-ep)*LP)./(h0+LP*thetap);

% horizontal velocity across the retrowedge
vretro=(hm*vm+(UR-er)*xr)./(hR+(LR-xr)*thetar);
 
% vertical velocity on the retrowedge surface
urs=er-(hm*vm+(UR-er)*xr)*alphar./((hR+(LR-xr)*thetar)); 
 
% elevation of retrowedge surface
zsr=elevm-xr*alphar;
 
% elevation of the base of the retrowedge
zbr=elevm-hm+xr*betar;

 %% integrate particle path

   
 np1=40; %number of particle paths for frontal accretion (was at 40)
 np2=10; %number of particle paths for underplating
 np=np1+np2;
 xp=zeros(np,1000);
 zp=zeros(np,1000);
 ntime=1000.;
 delt=10*(LP+LR)/(ntime*vp);
 
 % Horizontal velocity and Vertical velocity (for output into table)
 vHorz = zeros(np,ntime);
 vVert = zeros(np,ntime);
 
 % Depth between each frontally accreted particle
 delz=h0/(np1+1);
 
 % Horizontal distance between each underplated particle
 delx=(LP+LR)/(np2);
 
 % Starting x and z values for frontally accreted particles
 for ip=1:np1
  xp(ip,1)=0.0;
  zp(ip,1)=ip*delz;
 end
 
 % Starting x and z values for underplated particles
 for ip=np1+1:np
  xp(ip,1)=(ip-np1-1)*delx;
  if xp(ip,1) < LP
      zp(ip,1)=-xp(ip,1)*betap;
  else
      zp(ip,1)=-LP*betap+(xp(ip,1)-LP)*betar;
  end
  
 end
 
 for ip=1:np1+np2
   
    
    for it= 1:ntime-1
    
    % For particles within the Prowedge, perform the following loop
     if xp(ip,it) < LP
         % Equation for general horizontal velocity
         vpart=(h0*vp+(UP-ep)*xp(ip,it))/(h0+xp(ip,it)*thetap);
         vHorz(ip,it) = (h0*vp+(UP-ep)*xp(ip,it))/(h0+xp(ip,it)*thetap);
         
         % Equation for full vertical velocity that interpolates between the surface and base
         upart=UP-vpart*betap+(ep-UP+vpart*(alphap+betap))*(zp(ip,it)+xp(ip,it)*betap)/(h0+xp(ip,it)*(alphap+betap));
         vVert(ip,it) = UP-vpart*betap+(ep-UP+vpart*(alphap+betap))*(zp(ip,it)+xp(ip,it)*betap)/(h0+xp(ip,it)*(alphap+betap));
        
         xp(ip,it+1)=xp(ip,it)+delt*vpart;
         zp(ip,it+1)=zp(ip,it)+delt*upart;
         
         if zp(ip,it+1) > (h0+xp(ip,it+1)*alphap)
              zp(ip,it+1)=h0+xp(ip,it+1)*alphap;
         end
         
     % For particles within the Retrowedge, perform the following loop
     elseif xp(ip,it) < LP+LR
         % Equation for general horizontal velocity
         vpart=(hm*vm+(UR-er)*(xp(ip,it)-LP))/(hR+(LR-xp(ip,it)+LP)*thetar);
         vHorz(ip,it) = (hm*vm+(UR-er)*(xp(ip,it)-LP))/(hR+(LR-xp(ip,it)+LP)*thetar);
         
         % Equation for full vertical velocity that interpolates between the surface and base
         upart=UR+vpart*betar+(er-UR-vpart*(alphar+betar))*(zp(ip,it)+LP*betap-(xp(ip,it)-LP)*betar)/(h0+LP*thetap-(xp(ip,it)-LP)*thetar);
         vVert(ip,it) = UR+vpart*betar+(er-UR-vpart*(alphar+betar))*(zp(ip,it)+LP*betap-(xp(ip,it)-LP)*betar)/(h0+LP*thetap-(xp(ip,it)-LP)*thetar);
         
         xp(ip,it+1)=xp(ip,it)+delt*vpart;
         zp(ip,it+1)=zp(ip,it)+delt*upart;
         
         if zp(ip,it+1) > (h0+LP*alphap-(xp(ip,it+1)-LP)*alphar)
              zp(ip,it+1)=(h0+LP*alphap-(xp(ip,it+1)-LP)*alphar);
         end
     else
         xp(ip,it+1)=xp(ip,it);
         zp(ip,it+1)=zp(ip,it);
     end
     
    end
     
     
 end
 % make isoage dots
 freq=1; %frequency of isoage dots
 icount=0;
 xpspar=zeros(np,1000);
 zpspar=zeros(np,1000);

%xpspar and zpspar dictate where along material path the particle is plotted 
   for it=1:ntime
       m=mod(it,freq);
       if m == 0
           icount=icount+1;
           for ip=1:np
               xpspar(ip,icount)=xp(ip,it);
               zpspar(ip,icount)=zp(ip,it);
           end
       end
   end
           
          
%% Set variables for X and Y coordinates for closure depths along wedge

%Closure depth as a scalar
AHescalar = -2;
AFTscalar = -3.3;
ZHescalar = -6.;

%Closure depth that parallels topography
AHeCD = [zs+ AHescalar, zsr + AHescalar];
AFTCD = [zs+ AFTscalar, zsr + AFTscalar];
ZHeCD = [zs+ ZHescalar, zsr + ZHescalar];

%% Map distances onto actual x distances along wedge 

dist = LP+LR - xp;

% Depths of Material Paths
RockDepth = zp - h0;

%Clip Rock Depth when path reaches end of wedge
RockDepthClip = zeros(np,[]);

% Dummy variable for clipped distance
distClip = nan(np,ntime);

for ip = 1:np
    d1 = unique(dist(ip,:),'stable');
    d2 = unique(RockDepth(ip,:), 'stable');
    distClip(ip,1:length(d1)) = d1;
    distClip(isnan(distClip)) = [];
    RockDepthClip(ip,1:length(d2)) = d2;
end

% Common values between rows and within a single column indicate that samples are now traveling along
% the surface, so remove these surface points points 
for it = 1:length(distClip)
    for ip = np1:-1:1
        [r,~] = find(RockDepthClip(1:ip-1,it) == RockDepthClip(ip,it));
        if isempty(r)
            continue
        else  
            RockDepthClip(r,it) = nan;
            RockDepthClip(ip,it) = nan;
        end 
    end
end

%Manually set place at which first particle path hits the surface
RockDepthClip(np1,5:7) = nan;



%% Map closure depths onto distances along wedge

pro = LR+LP - x;
retro =LR-xr;

dist2 = [pro, retro];
dist2 = fliplr(unique(dist2));

AHedepth = zeros(np,ntime);
AFTdepth = zeros(np,ntime);
ZHedepth = zeros(np,ntime);

tol = 0.01;

% Create matrix of closure depths for all material paths
for ip = 1:np
    for p = 1:ntime       
          [~,col] =  find(abs((dist(ip, p)- dist2) <tol));
        if length(col) > 1  
           [~,col] = min(abs(dist(ip,p) - dist2(col)));
           
        end
        AHedepth(ip,p) = AHeCD(col);
        AFTdepth(ip,p) = AFTCD(col);
        ZHedepth(ip,p) = ZHeCD(col); 
        
    end
    
end

%% Calculate distance along wedge at which point sample passes above best-fit closure depth
% Note that calculations of CDy and CDx variables are only needed for plotting the data on the wedge

% Set tolerance for closure depth calculations
tol = 0.25;

% Record column number (x dist) where sample hits closure depth
AHecolInit = zeros(np,1);
AFTcolInit = zeros(np,1);
ZHecolInit = zeros(np,1);

%Dummy variable to catch best-fit closure depth from matrix for each material path and thermochronometer
AHeCDy = zeros(np,1);
AftCDy = zeros(np,1);
ZHeCDy = zeros(np,1);

% X value at which particle crosses best-fit closure depth
AHeCDx = zeros(np,1);
AFTCDx = zeros(np,1);
ZHeCDx = zeros(np,1);

% thermochronometer
for ip = 1:np
 
    %AHe
    column1 = find(abs(AHedepth(ip,:)-RockDepth(ip,:)) <tol);
    value1 = AHedepth(ip,:)-RockDepth(ip,:);   
    value1(value1<0)=nan;
    if isempty(column1) 
        AHeCDx(ip) = nan;  
        AHeCDy(ip) = nan;
        AHecolInit(ip) = nan;
    elseif isnan(value1)
        AHeCDx(ip) = nan;  
        AHeCDy(ip) = nan;
        AHecolInit(ip) = nan;
    elseif length(column1) > 1 
       [~,c] = min(value1);
       if RockDepth(ip,c) > RockDepth(ip,c+1)
            A2 = sort(value1);
            out = A2(2);
            [~,c] = find(value1 == out,1, 'first');
        else
        end
       AHeCDx(ip) = LP+LR - xp(ip,c);
       %AHeX(ip) = xp(ip,c);
       AHeCDy(ip) = AHedepth(ip,c);
       AHecolInit(ip) = c;
    else
       AHeCDx(ip) = LP+LR - xp(ip, column1);
       %AHeX(ip) = xp(ip, column1);
       AHeCDy(ip) = AHedepth(ip,column1);
       AHecolInit(ip) = column1;
    end
    
    %AFT
    column2 = find(abs(AFTdepth(ip,:)-RockDepth(ip,:)) <tol);
    value2 = AFTdepth(ip,:)-RockDepth(ip,:);
    value2(value2<0) =nan;
    if isempty(column2) 
        AFTCDx(ip) = nan;   
        AftCDy(ip) = nan;
        AFTcolInit(ip) = nan;
    elseif  isnan(value2)
        AFTCDx(ip) = nan;   
        AftCDy(ip) = nan;
        AFTcolInit(ip) = nan;
    elseif length(column2) > 1 
       [~,c] = min(value2);
       if RockDepth(ip,c) > RockDepth(ip,c+1)
            A2 = sort(value2);
            out = A2(2);
            [~,c] = find(value2 == out, 1, 'first');
        else
        end
       AFTCDx(ip) = LP+LR - xp(ip,c);
       %AFTX(ip) = xp(ip,c); 
       AftCDy(ip) = AFTdepth(ip,c);
       AFTcolInit(ip) = c;

    else
       AFTCDx(ip) = LP+LR - xp(ip, column2);
       %AFTX(ip) = xp(ip, column2);
       AftCDy(ip) = AFTdepth(ip,column2);
       AFTcolInit(ip) = column2;
    end
    
    %ZHe
    column3 = find(abs(ZHedepth(ip,:)-RockDepth(ip,:)) <tol);
    value3 = ZHedepth(ip,:)-RockDepth(ip,:); 
    value3(value3<0) =nan;
    if isempty(column3) 
        ZHeCDx(ip) = nan; 
        ZHeCDy(ip) = nan;
        ZHecolInit(ip) = nan;
    elseif isnan(value3)
        ZHeCDx(ip) = nan; 
        ZHeCDy(ip) = nan;
        ZHecolInit(ip) = nan;    
    elseif length(column3) > 1 
           [~,c] = min(value3);
           if RockDepth(ip,c) > RockDepth(ip,c+1)
            A2 = sort(value3);
            out = A2(2);
            [~,c] = find(value3 == out, 1, 'first');
           else
           end
        ZHeCDx(ip) = LP+LR - xp(ip,c);
       %ZHeX(ip) = xp(ip,c);
       ZHeCDy(ip) = ZHedepth(ip,c);
       ZHecolInit(ip) = c;
    else
       ZHeCDx(ip) = LP+LR - xp(ip,column3);
       %ZHeX(ip) = xp(ip,column3);
       ZHeCDy(ip) = ZHedepth(ip, column3);
       ZHecolInit(ip) = column3;
    end
end

%% Calculate X and Y Distances between CD and wedge surface

t = zeros(np, length(zpspar));
s = zeros(np, length(zpspar));

% Record column number (x dist) where sample hit surface
ColFinal = zeros(np,1);

% X Value at which rock reaches surface
SurfX = zeros(np,1);

xval = (xpspar(np1,:));
[~,col] = unique(xval, 'stable');
SurfZ = zpspar(np1,col);
SurfZ = SurfZ -h0;
SurfZ = reshape(SurfZ, length(distClip),1);

SurfZrock = zeros(np,1);

% maximum burial depth
Maxburial = zeros(np,1);
xPosition = zeros(np,1);
zPosition = zeros(np,1);

c = zeros(np,1);
   
for ip = 1:np
    for p = 1:length(SurfZ)
    [~,c(ip)] = find(~isnan(RockDepthClip(ip,:)), 1, 'last');
        if c(ip) <=length(distClip)
            SurfX(ip) = distClip(ip,c(ip));
            SurfZrock(ip) = SurfZ(c(ip));
            ColFinal(ip) = c(ip);
        elseif isempty(c(ip))
            [~,c(ip)] = find(~iszero(RockDepthClip(ip,:), 1, 'last'));
            if isempty(c(ip))
                ColFinal(ip) = nan;
            else
                SurfX(ip) = distClip(ip,c(ip));
                SurfZrock(ip) = SurfZ(c(ip));
                ColFinal(ip) = c(ip);   
            end
        elseif c(ip) >length(distClip)
            ColFinal(ip) = nan;
        end
    end
end

for ip = np1+1:np
    if zpspar(ip,c(ip))-h0 < SurfZ(ip)
        SurfX(ip) = nan;
        SurfZrock(ip) = nan;
    else
    end
end

for ip = 2:np
    if ColFinal(ip) == ColFinal(ip-1)
        ColFinal(ip-1) = nan;
    elseif isnan(ColFinal(np-1))
        ColFinal(np) = nan;
    else
    end

if isnan(SurfX(np-1))
    SurfX(np) = nan;
    SurfZrock(np) = nan;
else
end
end


%% Calculate time that particle travels as a function of its distance and velocity
AHeVel = zeros(np,ntime);
AFTVel = zeros(np,ntime);
ZHeVel = zeros(np,ntime);

AHedist = zeros(np,ntime);
AFTdist = zeros(np,ntime);
ZHedist = zeros(np,ntime);

AHetime = zeros(np,ntime);
AFTtime = zeros(np,ntime);
ZHetime = zeros(np,ntime);

AHetime(1:np,1) = nan;
AFTtime(1:np,1) = nan;
ZHetime(1:np,1) = nan;

AHetotTime = zeros(np,1);
AFTtotTime = zeros(np,1);
ZHetotTime = zeros(np,1);

AHetotDist = zeros(np,1);
AFTtotDist = zeros(np,1);
ZHetotDist = zeros(np,1);

% Calculate Time from distances and velocities 
for ip = 1:np
    for it = 2:ntime + 1
        if it-1 >= AHecolInit(ip) && it-1 <= ColFinal(ip)
            distx = xp(ip,it-1) - xp(ip,it);
            disty = zp(ip,it-1) - zp(ip,it);
            AHedist(ip,it-1) = sqrt(distx^2 + disty^2);
            AHetotDist(ip) = nansum(AHedist(ip,:));
            
            velHorz = vHorz(ip,it-1);
            velVert = vVert(ip,it-1);
            AHeVel(ip,it-1) = (sqrt(velVert^2 + velHorz^2));
            
            AHetime(ip,it-1) = (AHedist(ip,it-1)/AHeVel(ip,it-1));
        else
            AHedist(ip,it-1) = nan;
            AHeVel(ip,it-1) =nan;
            AHetime(ip,it-1) = nan;
        end
        AHetime(isinf(AHetime)) = nan;
            
        if it-1 >= AFTcolInit(ip) && it-1 <= ColFinal(ip)
            distx = xp(ip,it-1) - xp(ip,it);
            disty = zp(ip,it-1) - zp(ip,it);
            AFTdist(ip,it-1) = sqrt(distx^2 + disty^2);
            AFTtotDist(ip) = nansum(AFTdist(ip,:));
            
            velHorz = vHorz(ip,it-1);
            velVert = vVert(ip,it-1);
            AFTVel(ip,it-1) = sqrt(velVert^2 + velHorz^2);
            
            AFTtime(ip,it-1) = AFTdist(ip,it-1)/AFTVel(ip,it-1);
        else
            AFTdist(ip,it-1) = nan;
            AFTVel(ip,it-1) =nan;
            AFTtime(ip,it-1) = nan;
        end
         AFTtime(isinf(AFTtime)) = nan;
         
         if it-1 >= ZHecolInit(ip) && it-1 <= ColFinal(ip)   
            distx = xp(ip,it-1) - xp(ip,it);
            disty = zp(ip,it-1) - zp(ip,it);
            ZHedist(ip,it-1) = sqrt(distx^2 + disty^2);
            ZHetotDist(ip) = nansum(ZHedist(ip,:));
            
            velHorz = vHorz(ip,it-1);
            velVert = vVert(ip,it-1);
            ZHeVel(ip,it-1) = sqrt(velVert^2 + velHorz^2);
            
            ZHetime(ip,it-1) = ZHedist(ip,it-1)/ZHeVel(ip,it-1);
         else
             ZHedist(ip,it-1) = nan;
             ZHeVel(ip,it-1) = nan;
             ZHetime(ip,it-1) = nan;
         end
        ZHetime(isinf(ZHetime)) = nan;
    end
end

% Flip time variable so that it matches pattern of X distances
    AHetime = fliplr(AHetime);
    AFTtime = fliplr(AFTtime);
    ZHetime = fliplr(ZHetime);
    
    erateAHe = zeros(np,1);
    erateAFT = zeros(np,1);
    erateZHe = zeros(np,1);

for ip = 1:np
    if any(~isnan(AHetime(ip,:)))
        AHetotTime(ip) = nansum(AHetime(ip,:),2);
        
    else
         AHetotTime(ip) = nan; 
    end
    
    if any(~isnan(AFTtime(ip,:)))
        AFTtotTime(ip) = nansum(AFTtime(ip,:),2);
            
    else
        AFTtotTime(ip) = nan; 
    end
    
    if any(~isnan(ZHetime(ip,:)))
         ZHetotTime(ip) = nansum(ZHetime(ip,:),2);    
    else
        ZHetotTime(ip) = nan;
    end

  % Calculate an erosion rate for each thermochronometer from cooling age and distance to surface    
    erateAHe(ip) = abs(AHescalar)/AHetotTime(ip);
    erateAFT(ip) = abs(AFTscalar)/AFTtotTime(ip);
    erateZHe(ip) = abs(ZHescalar)/ZHetotTime(ip);
    
    Maxburial(ip) = max(zp(ip));
    xPosition(ip) = SurfX(ip);
    zPosition(ip) = SurfZrock(ip);
end


%% Horizontal Velocities, Uplift, and Wedge Material Path Plots

figure(1)
set(gcf, 'Position',  [600, 600,1000, 900])
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
 subplot(3,1,1)
 plot(LP+LR-x,v, LR-xr,vretro, 'Linewidth', 2);
 %xlabel('Distance [km]'); 
 ylabel('Horizontal Velocity (mm/yr)')
 axis([0 LP+LR  0 20])
 set(gca,'FontSize',14)
 set(gca,'XTickLabel',{' '})
 
 subplot(3,1,2)
 xl(1)=0;
 xl(2)=LP+LR;
 yl(1)=0;
 yl(2)=0;
 plot(LP+LR-x,us,LR-xr,urs, xl,yl,'k', 'Linewidth', 2)
 %xlabel('Distance (km)'); 
 ylabel('Uplift Rate (mm/yr)')
 axis([0 LP+LR  -2 2.])
 set(gca,'FontSize',14)
 set(gca,'XTickLabel',{' '})
 
subplot(3,1,3)
 
 xlabel('Distance (km)'); ylabel('Elevation (km)')
 axis([0 LP+LR  -60 5])
 set(gca,'FontSize',14)

 % comment out for third plot
hold
 for ip=1:np
    t(ip,:) =  zpspar(ip,:)-h0;   
    s(ip,:) = LP+LR-xpspar(ip,:);
     for n = 1:20:ntime
         if ip < np1 && mod(ip,5) == 0
             plot(LP+LR-xp(ip,:),zp(ip,:)-h0, 'LineWidth', 1)
             scatter(s(ip,n),t(ip,n),'filled')
         elseif ip > np1 
             plot(LP+LR-xp(ip,:),zp(ip,:)-h0, 'LineWidth', 1)
             scatter(s(ip,n),t(ip,n),'filled')
         else 
         end
     end
 end
 plot(LP+LR-x,zs,'k-',LR-xr,zsr,'k-',LP+LR-x,zb,'k-',LR-xr,zbr,'k-', 'LineWidth', 3)
 
 print('-painters', '-dpdf','Wedge_Material_Paths')

%% Predicted Age vs Distance Plot
figure(2)

set(gcf, 'Position',  [600, 600,1000, 300])
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));


scatter(SurfX,AHetotTime,125,'^', 'Markerfacecolor', [0 0 1],'Markeredgecolor', 'k');hold on
row = find(~isnan(AHetotTime), 1, 'first');
text(max(SurfX), max(AHetotTime),'AHe', 'FontSize', 14)

scatter(SurfX,AFTtotTime,125,'d','Markerfacecolor', [1 0 0], 'Markeredgecolor', 'k');
row = find(~isnan(AFTtotTime), 1, 'first');
text(max(SurfX),max(AFTtotTime) ,'AFT', 'FontSize', 14)

scatter(SurfX,ZHetotTime, 125,'o','Markerfacecolor', [0 1 0],  'Markeredgecolor', 'k')
row = find(~isnan(ZHetotTime), 1, 'first');
text(max(SurfX), max(ZHetotTime),'ZHe', 'FontSize', 14)

xlabel('Distance (km)', 'FontSize', 14); ylabel('Predicted Age Ma)', 'FontSize', 14)
axis([0 LP+LR  0 15])
set(gca,'FontSize',14)

%% Reset Thermochronometers Plot
figure(3)
set(gcf, 'Position',  [100, 100,1000, 600])

plot(LP+LR-x,zs,'k-',LR-xr,zsr,'k-',LP+LR-x,zb,'k-',LR-xr,zbr,'k-', 'LineWidth', 3); hold on

plot(LP+LR-x,zs+ AHescalar,'--k',LR-xr,zsr+ AHescalar,'--k','LineWidth', 1.5); 
text(101, -1.8,'AHe', 'FontSize', 14)

plot(LP+LR-x,zs+ AFTscalar,'--k',LR-xr,zsr+ AFTscalar,'--k','LineWidth', 1.5); 
text(101, -3.1,'AFT', 'FontSize', 14)

plot(LP+LR-x,zs+ ZHescalar,'--k',LR-xr,zsr+ ZHescalar,'--k','LineWidth', 1.5); 
text(101, -5.8,'ZHe', 'FontSize', 14)

xlabel('Distance (km)', 'FontSize', 14); ylabel('Elevation (km)', 'FontSize', 14)
axis([0 LP+LR  -8 2.5])
set(gca,'FontSize',14)

%Extract closure temperatures for all thermochrons for color mapping
totTime = [AHetotTime; AFTtotTime; ZHetotTime];

 for ip=1:np 
 %AHe
    if ~isnan(AHetotTime(ip)) && SurfX(ip) >0 && AHeCDx(ip) < 40
        plot(LP+LR-xp(ip,:),zp(ip,:)-h0, '-k', 'Linewidth', 1); hold on
        scatter(AHeCDx(ip),AHeCDy(ip),125,'^', 'Markerfacecolor', [0.75 0.75 0.75],'Markeredgecolor', 'k');
        str = num2str(AHetotTime(ip), '%.1f');
%        text(AHeCDx(ip), AHeCDy(ip)+0.4, str)
    elseif ~isnan(AHetotTime(ip)) &&  SurfX(ip) >0 && AHeCDx(ip) > 40
        plot(LP+LR-xp(ip,:),zp(ip,:)-h0, '-k', 'Linewidth', 1); 
        scatter(AHeCDx(ip),AHeCDy(ip),125,'^', 'Markerfacecolor', [0.6 0.6 0.6],'Markeredgecolor', 'k'); hold on
        str = num2str(AHetotTime(ip), '%.1f');
%        text(AHeCDx(ip), AHeCDy(ip)+0.4, str)
    else
        plot(LP+LR-xp(ip,:),zp(ip,:)-h0, 'Color', [0.5 0.5 0.5])
    end
  
 % AFT
    if ~isnan(AFTtotTime(ip)) && SurfX(ip) >0 && AFTCDx(ip) < 40
        scatter(AFTCDx(ip),AftCDy(ip),125,'d','Markerfacecolor', [0.4 0.4 0.4], 'Markeredgecolor', 'k');
        str = num2str(AFTtotTime(ip),'%.1f');
    %    text(AFTCDx(ip), AftCDy(ip)+0.4, str)
    elseif ~isnan(AFTtotTime(ip)) && ~isinf(AFTtotTime(ip)) && AFTCDx(ip) > 40
        scatter(AFTCDx(ip),AftCDy(ip),125,'d','Markerfacecolor', [0.4 0.4 0.4], 'Markeredgecolor', 'k');
        str = num2str(AFTtotTime(ip),'%.1f');
    %    text(AFTCDx(ip), AftCDy(ip)+0.4, str)
    else
     
    end
 
%ZHe
    if ~isnan(ZHetotTime(ip)) && SurfX(ip) >0 && ZHeCDx(ip) <40
        scatter(ZHeCDx(ip),ZHeCDy(ip), 125,'o','Markerfacecolor', [0.2 0.2 0.2],  'Markeredgecolor', 'k')
        str = num2str(ZHetotTime(ip),'%.1f');
    %    text(ZHeCDx(ip), ZHeCDy(ip)+0.4, str, 'FontSize', 10)
    elseif ~isnan(ZHetotTime(ip)) && ~isinf(ZHetotTime(ip)) && ZHeCDx(ip) >40
        scatter(ZHeCDx(ip),ZHeCDy(ip), 125, 'o','Markerfacecolor', [0.2 0.2 0.2], 'Markeredgecolor', 'k')
        str = num2str(ZHetotTime(ip),'%.1f');
    %   text(ZHeCDx(ip), ZHeCDy(ip)+0.4, str, 'FontSize', 10)
    else
    end
 
scatter(SurfX(ip),SurfZrock(ip),100,'s','Markerfacecolor', [0.8 0.8 0.8], 'Markeredgecolor', 'k');

 end

for i = 1:length(Maxburial)
    if isnan(zPosition(i)) || zPosition(i) < 0
        Maxburial(i) = NaN;
        xPosition(i) = NaN;
    else
    end
end
    

print('-painters','Wedge_Cooling_Ages_w_Labels', '-dpdf')


