%% Wedge Kinematics

%Given a wedge with a frontal flux in, underplating, erosion and rearward flux out. Assume steady
%state, so that the fluxes balance.

%Change parameters in Lines 26-35.

%The code outputs 4 plots: 

%(1)Wedge material paths for the full orogen with predicted uplift rates 
%   and horizontal velocities across the wedge
%(2)The predicted ages of each thermochronometer across the orogenic wedge
%(3)Material paths for reset and unreset thermochronometers in the upper 8 
%km of the orogenic wedge
%(4)Maximum burial depths for reset and unreset particles

clear all
close all

cd 'C:\Users\erlanger\Desktop\Matlab\Wedge_Model'

addpath 'C:\Users\erlanger\Desktop\Matlab\Wedge_Model'
addpath 'C:\Users\erlanger\Desktop\Matlab'

%%  Set the following parameters:
h0=10.;  % incoming section for frontal accretion In %km
h1=10.;  % incoming section for prowedge basal accretion
h2=0.;  % incoming section for retrowedge basal accretion
ep=0.6;  % pro-wedge erosion rate%km/Ma 
er=0.3;  % retro-wedge erosion rate 
LP=60.; % length of Pro-wedge % km
LR=40.;  % Length of retro-Wedge
vp=10;  %Slab retreat rate %km/Ma  (Normally signifies convergence rate but Apennines has negligible convergence)
hm=56.; % Max Crustal Thickness
elevm=2;  % Maximum Elevation (km)

%% Wedge geometry

x = linspace(0,LP,1001);
u = zeros(1001,101);
xr = linspace(0,LR,1001);


%Prowedge Slopes
alphap = elevm/LP;
betap =(hm-elevm-h0)/LP;
thetap = alphap+betap;

% elevation of prowedge surface 
zs = x*alphap; 

% elevation of the base of the prowedge
zb = -h0-x*betap;

%determine thickness of outgoing crust hR 
hR = h0+h1+h2-ep*LP/vp-er*LR/vp;
retrocrust = hR;

%Retro-wedge Slopes
alphar =(elevm)/LR;
betar =(hm-elevm-hR)/LR;
thetar = alphar+betar;

%calculate maximum crustal thickness from prowedge slopes
hm = h0+LP*(alphap+betap);
elevm = LP*(alphap);

z = linspace(0,hm,101);

% elevation of retrowedge surface
zsr = elevm-xr*alphar;
 
% elevation of the base of the retrowedge
zbr=elevm-hm+xr*betar;


%% Prowedge Field Velocities

% Underplating velocity on prowedge
UP = vp*h1/LP;

%find horizontal velocity across pro-wedge
vpro =(h0*vp+(UP-ep)*x)./(h0+x*thetap);

%find vertical velocity on the prowedge surface (defined as the erosion velocity plus the motion to carry the point
%along the surface)
us=ep+((h0*vp+x*(UP-ep))*alphap)./((h0+x*thetap));


%% Retrowedge Field Velocities

%Underplating velocity on retrowedge
UR = vp*h2/LR;

% velocity at the divide
vm =(h0*vp+(UP-ep)*LP)./(h0+LP*thetap);

% horizontal velocity across the retrowedge
vretro = (hm*vm+(UR-er)*xr)./(hR+(LR-xr)*thetar);
 
% vertical velocity on the retrowedge surface
urs = er-(hm*vm+(UR-er)*xr)*alphar./((hR+(LR-xr)*thetar)); 
 

 %% Integrate Particle Path

 np1 = 35; %number of particle paths for frontal accretion (was at 40)
 np2 = 10; %number of particle paths for underplating
 np = np1+np2;
 xp = zeros(np,1000);
 zp = zeros(np,1000);
 ntime = 1000.;  
 delt = 10*(LP+LR)/(ntime*vp);
 
 % Horizontal velocity and Vertical velocity (for output into table)
 vHorz = zeros(np,ntime);
 vVert = zeros(np,ntime);
 
 % Depth between each frontally accreted particle
 delz = h0/(np1+1);
 
 % Horizontal distance between each underplated particle
 delx = (LP+LR)/(np2);
 
 % Starting x and z values for frontally accreted particles
 for ip = 1:np1
  xp(ip,1) = 0.0;
  zp(ip,1) = ip*delz;
 end
 
 % Starting x and z values for underplated particles
 for ip = np1+1:np
  xp(ip,1)=(ip-np1-1)*delx;
  if xp(ip,1) < LP
      zp(ip,1)=-xp(ip,1)*betap;
  else
      zp(ip,1)=-LP*betap+(xp(ip,1)-LP)*betar;
  end 
 end
 
 for ip = 1:np  
    for it = 1:ntime-1
    
    % For particles within the Prowedge
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
         
     % For particles within the Retrowedge
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
 freq = 1; %frequency of isoage dots
 icount = 0;
 xpspar = zeros(np,1000);  
 zpspar = zeros(np,1000);

%xpspar and zpspar represent wedge coordinates along material path where 
%the particle is plotted
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
                   
%% Map distances onto actual x distances along wedge 

% Depths of Material Paths
rockDepth = zpspar - h0;

%Clip Rock Depth when path reaches end of wedge
rockDepthClip = zeros(np,[]);

% Material path (x direction) along wedge
xPath = LP+LR - xpspar;

% clipped X Distance along wedge
distClip = nan(np,ntime);

for ip = 1:np
    d1 = unique(xPath(ip,:),'stable');
    d1(isnan(d1)) = [];
    
    depth = unique(rockDepth(ip,:), 'stable');
    depth(isnan(depth))=[];
    
    distClip(ip,1:length(d1)) = d1;
    distClip(isnan(distClip)) = [];
    
    rockDepthClip(ip,1:length(depth)) = depth;
end

%% Find where particle path hits the surface of the wedge

% Matrix column that corresponds to intersection of the particle with the wedge surface 
colFinal =  zeros (np,1);

% Depth at which particle intersects the surface
surfZrock = zeros(np,1);

% Distance along wedge at which rock reaches surface
SurfX = zeros(np,1);

%Define (x,z) coordinates of wedge surface
prowedge = [x;zs]';
retrowedge = [fliplr(100-xr);zsr]';

%Define (x,z) particle coordinates
particles = cell(np,1);
for ip = 1:np
    particles{ip} = [xpspar(ip,:);zpspar(ip,:)-h0]'; 
end

%Find intersection of wedge surface with particle coordinates
%Prowedge
for ip = 1:np
    flag = 0;
    for it = 1:length(prowedge)
        [r,~] = find(round(particles{ip,1}(:,2),2) == round(prowedge(it,2),2), 1, 'first');
        if isempty(r)
               SurfX(ip) = nan;
               surfZrock(ip) = nan;
            continue
        elseif abs(round(particles{ip,1}(r,1),2) - round(prowedge(it,1),2)) < 0.5
            SurfX(ip) = 100-(particles{ip,1}(r,1));
            surfZrock(ip) = particles{ip,1}(r,2);
            colFinal(ip) = r;
            flag = 1;    
            break
        end
        if(flag==1)
            break
        end
    end
end

%Retrowedge
for ip = 1:np
    flag = 0;
    if ~isnan(SurfX(ip)) 
        continue
    else
        for it = 1:length(retrowedge) 
            [r,~] = find(round(particles{ip,1}(:,2),2) == round(retrowedge(it,2),2));  
            if ~isempty(r)   
                for i = 1: length(r)
                    if abs(round(particles{ip,1}(r(i),1),2) - round(retrowedge(it,1),2)) < 0.02 
                        SurfX(ip) = 100-(particles{ip,1}(r(i),1));
                        surfZrock(ip) = particles{ip,1}(r(i),2);
                        colFinal(ip) = r(i);
                        flag = 1;  
                        break
                    else
                        continue
                    end
                end
            else
            end
            if(flag==1)
                break
            end
        end
    end
end

%% Closure Depth

%Closure depth as a scalar calculated from closure temperatures (see
%Farley, (2000) for AHe and ZHe and Wagner and Van den Haute (1992) for
%AFT) and an average geothermal gradient (from heat flow measurements) 
% = 36.4 degree C/km
AHescalar = -1.9;
AFTscalar = -3.0;
ZHescalar = -4.9;

% Closure depth that parallels wedge topography
AHeCDwedge = [zs+ AHescalar, zsr + AHescalar];
AFTCDwedge = [zs+ AFTscalar, zsr + AFTscalar];
ZHeCDwedge = [zs+ ZHescalar, zsr + ZHescalar];


%% Map closure depths onto distances along wedge

%Prowedge
pro = LR+LP-x;

%Retrowedge
retro =LR-xr;

%Full wedge
wedgeDist = [pro, retro];
wedgeDist = fliplr(unique(wedgeDist));

AHeCD = zeros(np,ntime);
AFTCD = zeros(np,ntime);
ZHeCD = zeros(np,ntime);

%Set tolerance
tol = 0.01;

% Calculate Closure depths for all particles
for ip = 1:np
    for p = 1:ntime  
        % Find column within xPath row that most closely matches value in
        % wedgeDist array
        [~,col1] =  find((xPath(ip, p)- wedgeDist) <tol);
        if length(col1) > 1  
           [~,col1] = min(abs(xPath(ip,p) - wedgeDist(col1)));
        elseif isempty(col1)
            continue
        else
        end
        % Call column to produce matrix of correct closure depths that parallel wedge
        % topography
        AHeCD(ip,p) = AHeCDwedge(col1);
        AFTCD(ip,p) = AFTCDwedge(col1);
        ZHeCD(ip,p) = ZHeCDwedge(col1);    
    end       
end

%% Calculate distance along wedge at which point sample passes above closure depth
% Note that calculations of CDy and CDx variables are only needed for plotting the data on the wedge

% Set tolerance for closure depth calculations
tol = 0.25;

% Record column number (x dist) where sample hits closure depth
AHecol = zeros(np,1);
AFTcol = zeros(np,1);
ZHecol = zeros(np,1);

% Particle Closure Depth 
AHeCDy = zeros(np,1);
AftCDy = zeros(np,1);
ZHeCDy = zeros(np,1);

% X distance at which particle crosses closure depth
AHeCDx = zeros(np,1);
AFTCDx = zeros(np,1);
ZHeCDx = zeros(np,1);

% thermochronometer
for ip = 1:np
 
    %AHe
     col1 = find(abs(AHeCD(ip,:)-rockDepth(ip,:)) <tol);
     val1 = AHeCD(ip,:)-rockDepth(ip,:);   
     val1(val1<0)=nan;
    if isempty(col1)  
        AHeCDx(ip) = nan;  
        AHeCDy(ip) = nan;
        AHecol(ip) = nan;
    elseif isnan(val1)
        AHeCDx(ip) = nan;  
        AHeCDy(ip) = nan;
        AHecol(ip) = nan;
    elseif length(col1) > 1 
       [~,c] = min(val1);
       if rockDepth(ip,c) > rockDepth(ip,c+1)
            A2 = sort(val1);
            out = A2(2);
            [~,c] = find(val1 == out,1, 'first');
        else
        end
       AHeCDx(ip) = xPath(ip,c);
       AHeCDy(ip) = AHeCD(ip,c);
       AHecol(ip) = c;
    else
       AHeCDx(ip) = xPath(ip, col1);
       AHeCDy(ip) = AHeCD(ip,col1);
       AHecol(ip) = col1;
    end
    
    %AFT
    col2 = find(abs(AFTCD(ip,:)-rockDepth(ip,:)) <tol);
    val2 = AFTCD(ip,:)-rockDepth(ip,:);
    val2(val2<0) =nan;
    if isempty(col2)  
        AFTCDx(ip) = nan;   
        AftCDy(ip) = nan;
        AFTcol(ip) = nan;
    elseif isnan(val2)
        AFTCDx(ip) = nan;   
        AftCDy(ip) = nan;
        AFTcol(ip) = nan;
    elseif length(col2) > 1 
       [~,c] = min(val2);
       if rockDepth(ip,c) > rockDepth(ip,c+1)
            A2 = sort(val2);
            out = A2(2);
            [~,c] = find(val2 == out, 1, 'first');
        else
        end
       AFTCDx(ip) = xPath(ip,c);
       AftCDy(ip) = AFTCD(ip,c);
       AFTcol(ip) = c;
    else
       AFTCDx(ip) = xPath(ip, col2);
       AftCDy(ip) = AFTCD(ip,col2);
       AFTcol(ip) = col2;
    end
    
    %ZHe
    col3 = find(abs(ZHeCD(ip,:)-rockDepth(ip,:)) <tol);
    val3 = ZHeCD(ip,:)-rockDepth(ip,:); 
    val3(val3<0) =nan;
    if isempty(col3)  
        ZHeCDx(ip) = nan; 
        ZHeCDy(ip) = nan;
        ZHecol(ip) = nan;  
    elseif isnan(val3)
        ZHeCDx(ip) = nan; 
        ZHeCDy(ip) = nan;
        ZHecol(ip) = nan; 
    elseif length(col3) > 1 
        [~,c] = min(val3);
        if rockDepth(ip,c) > rockDepth(ip,c+1)
            A2 = sort(val3);
            out = A2(2);
            [~,c] = find(val3 == out, 1, 'first');
       else
       end
        ZHeCDx(ip) = xPath(ip,c);
       ZHeCDy(ip) = ZHeCD(ip,c);
       ZHecol(ip) = c;
    else
       ZHeCDx(ip) = xPath(ip,col3);
       ZHeCDy(ip) = ZHeCD(ip, col3);
       ZHecol(ip) = col3;
    end
end

%% Calculate time that particle travels as a function of its distance and velocity

%Particle Velocities
AHeVel = zeros(np,ntime);
AFTVel = zeros(np,ntime);
ZHeVel = zeros(np,ntime);

%Particle Distance
AHedist = zeros(np,ntime);
AFTdist = zeros(np,ntime);
ZHedist = zeros(np,ntime);

%Particle Time
AHetime = zeros(np,ntime);
AFTtime = zeros(np,ntime);
ZHetime = zeros(np,ntime);

%Particle Total Time
AHetotTime = zeros(np,1);
AFTtotTime = zeros(np,1);
ZHetotTime = zeros(np,1);

% Particle Total Distance
AHetotDist = zeros(np,1);
AFTtotDist = zeros(np,1);
ZHetotDist = zeros(np,1);

%Calculate Time from distances and velocities 
for ip = 1:np
    for it = 2:ntime + 1
        if it-1 >= AHecol(ip) && it-1 <= colFinal(ip)
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
            
        if it-1 >= AFTcol(ip) && it-1 <= colFinal(ip)
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
         
         if it-1 >= ZHecol(ip) && it-1 <= colFinal(ip)   
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

%Flip time variable so that it matches pattern of X distances
AHetime = fliplr(AHetime);
AFTtime = fliplr(AFTtime);
ZHetime = fliplr(ZHetime);

%Maximum burial depth
Maxburial = zeros(np,1);
xPosition = SurfX;
zPosition = surfZrock;
    
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
    
    %Maxburial is the most negative depth value (min value)
    Maxburial(ip) = min(rockDepth(ip,:));
    
    if isnan(xPosition(ip)) || xPosition(ip) <0
        Maxburial(ip) = nan;
    else
    end
end

Maxburial = abs(Maxburial);

%If Maxburial of particle is less than AHe closure depth or closure depth
%is nan, then SurfX is nan
for i = 1:np
    if Maxburial(i) < AHeCDy(i)
        SurfX(i) = nan;
    elseif isnan(AHeCDy(i))
        SurfX(i) = nan;
    end
end

%% Plot Data

t = zeros(np, length(zpspar));
s = zeros(np, length(zpspar));

figure(1)
set(gcf, 'Position',  [600, 600,1000, 900])
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
 
% Horizontal Velocity with Distance along the Wedge
subplot(3,1,1)
plot(pro,vpro, retro,vretro, 'Linewidth', 2);
ylabel('Horizontal Velocity (mm/yr)')
axis([0 LP+LR  0 20])
set(gca,'FontSize',14)
set(gca,'XTickLabel',{' '})

%Uplift Rate with Distance along the Wedge
subplot(3,1,2)
xl(1)=0;
xl(2)=LP+LR;
yl(1)=0;
yl(2)=0;
plot(pro,us,retro,urs, xl,yl,'k', 'Linewidth', 2)
ylabel('Uplift Rate (mm/yr)')
axis([0 LP+LR  -2 2.])
set(gca,'FontSize',14)
set(gca,'XTickLabel',{' '})
 
% Material motion paths through the wedge
subplot(3,1,3)
 
xlabel('Distance (km)'); ylabel('Elevation (km)')
axis([0 LP+LR  -60 5])
set(gca,'FontSize',14)

hold
 for ip=1:np
    t(ip,:) = zpspar(ip,:)-h0;   
    s(ip,:) = LP+LR-xpspar(ip,:);
     for n = 1:20:ntime
         if ip < np1 && mod(ip,5) == 0
             plot(xPath(ip,:),zp(ip,:)-h0, 'LineWidth', 1)
             scatter(s(ip,n),t(ip,n),'filled')
         elseif ip > np1 
             plot(xPath(ip,:),zp(ip,:)-h0, 'LineWidth', 1)
             scatter(s(ip,n),t(ip,n),'filled')
         else 
         end
     end
 end
 plot(pro,zs,'k-',retro,zsr,'k-',pro,zb,'k-',retro,zbr,'k-', 'LineWidth', 3)
 
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
text(max(SurfX), max(AHetotTime),'AHe', 'FontSize', 14)

scatter(SurfX,AFTtotTime,125,'d','Markerfacecolor', [1 0 0], 'Markeredgecolor', 'k');
text(max(SurfX),max(AFTtotTime) ,'AFT', 'FontSize', 14)

scatter(SurfX,ZHetotTime, 125,'o','Markerfacecolor', [0 1 0],  'Markeredgecolor', 'k')
text(max(SurfX), max(ZHetotTime),'ZHe', 'FontSize', 14)

xlabel('Distance (km)', 'FontSize', 14); ylabel('Predicted Age Ma)', 'FontSize', 14)
axis([0 LP+LR  0 15])
set(gca,'FontSize',14)

print('-painters', '-dpdf','Predicted_Age_vs_Distance')

%% Reset Thermochronometers Plot
figure(3)

set(gcf, 'Position',  [100, 100,1000, 600])

plot(pro,zs,'k-',retro,zsr,'k-',pro,zb,'k-',retro,zbr,'k-', 'LineWidth', 3); hold on

plot(pro,zs+ AHescalar,'-b',retro,zsr+ AHescalar,'-b','LineWidth', 1.5); 
text(101, -1.8,'AHe', 'FontSize', 14)

plot(pro,zs+ AFTscalar,'-r',retro,zsr+ AFTscalar,'-r','LineWidth', 1.5); 
text(101, -3.1,'AFT', 'FontSize', 14)

plot(pro,zs+ ZHescalar,'-g',retro,zsr+ ZHescalar,'-g','LineWidth', 1.5); 
text(101, -5.8,'ZHe', 'FontSize', 14)

xlabel('Distance (km)', 'FontSize', 14); ylabel('Elevation (km)', 'FontSize', 14)
axis([0 LP+LR  -8 2.5])
set(gca,'FontSize',14)

 for ip=1:np 
 %AHe
if ~isnan(AHetotTime(ip)) && SurfX(ip) >0 && AHeCDx(ip) < 40
    plot(xPath(ip,:),zp(ip,:)-h0, '-k', 'Linewidth', 1); hold on
    scatter(AHeCDx(ip),AHeCDy(ip),125,'^', 'Markerfacecolor', 'blue','Markeredgecolor', 'k');
    str = num2str(AHetotTime(ip), '%.1f');
    text(AHeCDx(ip), AHeCDy(ip)+0.4, str)
 elseif ~isnan(AHetotTime(ip)) &&  SurfX(ip) >0 && AHeCDx(ip) > 40
    plot(xPath(ip,:),zp(ip,:)-h0, '-k', 'Linewidth', 1); 
    scatter(AHeCDx(ip),AHeCDy(ip),125,'^', 'Markerfacecolor', 'blue','Markeredgecolor', 'k'); hold on
    str = num2str(AHetotTime(ip), '%.1f');
    text(AHeCDx(ip), AHeCDy(ip)+0.4, str)
else
    plot(xPath(ip,:),zp(ip,:)-h0, '--', 'Color', [0.5 0.5 0.5])
 end
  
 % AFT
 if ~isnan(AFTtotTime(ip)) && SurfX(ip) >0 && AFTCDx(ip) < 40
    scatter(AFTCDx(ip),AftCDy(ip),125,'d','Markerfacecolor', 'red', 'Markeredgecolor', 'k');
    str = num2str(AFTtotTime(ip),'%.1f');
    text(AFTCDx(ip), AftCDy(ip)+0.4, str)
 elseif ~isnan(AFTtotTime(ip)) && ~isinf(AFTtotTime(ip)) && AFTCDx(ip) > 40
    scatter(AFTCDx(ip),AftCDy(ip),125,'d','Markerfacecolor', 'red', 'Markeredgecolor', 'k');
    str = num2str(AFTtotTime(ip),'%.1f');
    text(AFTCDx(ip), AftCDy(ip)+0.4, str)
 else    
 end
 
%ZHe
 if ~isnan(ZHetotTime(ip)) && SurfX(ip) >0 && ZHeCDx(ip) <40
    scatter(ZHeCDx(ip),ZHeCDy(ip), 125,'o','Markerfacecolor', 'green',  'Markeredgecolor', 'k')
    str = num2str(ZHetotTime(ip),'%.1f');
    text(ZHeCDx(ip), ZHeCDy(ip)+0.4, str, 'FontSize', 10)
 elseif ~isnan(ZHetotTime(ip)) && ~isinf(ZHetotTime(ip)) && ZHeCDx(ip) >40
    scatter(ZHeCDx(ip),ZHeCDy(ip), 125, 'o','Markerfacecolor', 'green', 'Markeredgecolor', 'k')
    str = num2str(ZHetotTime(ip),'%.1f');
    text(ZHeCDx(ip), ZHeCDy(ip)+0.4, str, 'FontSize', 10)
else
 end
 
scatter(SurfX(ip),surfZrock(ip),100,'s','Markerfacecolor', [0.8 0.8 0.8], 'Markeredgecolor', 'k');

 end
 
print('-painters', '-dpdf','Reset_Thermochronometers', '-bestfit')

%% Predicted Maximum Burial Depth vs Distance Plot

figure(4)

set(gcf, 'Position',  [600, 600,1000, 300])
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));

for ip = 1:np
    if ~isnan(AHetotTime(ip))
        scatter(xPosition(ip),Maxburial(ip),'o','Markerfacecolor', 'k', 'Markeredgecolor', 'k'); hold on
    else
        scatter(xPosition(ip),Maxburial(ip),'o', 'Markeredgecolor', 'k')        
    end
end

xlabel('Distance (km)', 'FontSize', 14); ylabel('Maximum burial depth (km)', 'FontSize', 14)
axis([0 LP+LR 0 20])
set(gca,'FontSize',14)

print('-painters', '-dpdf','Predicted_Maximum_Burial')
