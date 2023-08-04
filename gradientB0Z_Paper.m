% Code to generate Figures 1 and 2 in Paper
% "The Effect of and Correction for Through-Slice Dephasing on 2D Gradient Echo Double Angle B1+ Mapping"

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2023

%% Phase Roll as a function of FA - Figure 1 (c)
clear;clc;close all
B1Bias = 1;
T1 = 900; %ms
T2 = 30; %ms
TE = 1:15; %ms %time from the centre of the RF pulse
kFactor = 2; %k*alpha=2alpha
FatSupress = 'FS';

% FA 5 until 160 degrees
sliceposition = -0.5:0.01:0.5; %cm
OffResSliceCentre = 0;%Hz
df_B0Zero = OffResSliceCentre*ones(1,length(sliceposition));
DeltaPhase = zeros(1,length(5:5:160));
cnt=0;
for inominalAngle = 5:5:160
    cnt=cnt+1;
[ratio_2alpha_alpha_B0Zero_VaryTE_FA65,~,~,mxy_alpha_B0Zero_VaryTE_FA65,mxy_2alpha_B0Zero_VaryTE_FA65,~,mxy_alpha_B0Zero_t0_FA65,mxy_2alpha_B0Zero_t0_FA65]= Mxy_te_Grefz(B1Bias,deg2rad(inominalAngle),kFactor,df_B0Zero,T1,T2,TE,sliceposition,FatSupress);
DeltaPhase(1,cnt) = abs(max(angle(mxy_alpha_B0Zero_t0_FA65(:,1))))+abs(min(angle(mxy_alpha_B0Zero_t0_FA65(:,1))));
end

fig=figure();
fig.Color = [1 1 1];
FA_Array = 5:5:160;
plot(FA_Array,DeltaPhase,'-','LineWidth',2);
xlabel('FA (degrees)')
ylabel('Phase Roll Extent (radians)');
grid on; 
xt = -0.15;
yt = 1.04;
str = {'(c)'};
text(xt,yt,str,'Units','normalized','fontsize',18,'FontName','Arial','FontWeight','bold')
set(gca,'FontSize',18,'FontName','Arial')

set(gca,'FontSize',18,'FontName','Arial')
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');
figName = 'PhaseRoll_FA5160';
%% MRM Letter Simulations: B0=0Hz; Constant B0, Linear B0gradientZ, Non-linear B0gradientZ

clc; clear;close all
%Constant Parameters 
T1 = 900; %ms
T2 = 30; %ms
TE = 1:15; %ms %time from the centre of the RF pulse
kFactor = 2; %k*alpha=2alpha

B1Bias = 1;
% sliceposition = -1:0.05:1; %cm
sliceposition = -1:0.01:1; %cm
nZPositions = length(sliceposition);
SliceThicknessFactor = 1; %it's 1 for 8mm slice thickness, 2 for 4 mm slice thickness
FatSupress = 'FS';
opts = optimset('Display','off');
B1Guess = 1;
shift = 0;
%======1. No Off-resonance: B0=0Hz========
OffResSliceCentre = 0;%Hz
% Small FA 5 degrees
nominalAngle = 5;
df_B0Zero = OffResSliceCentre*ones(1,length(sliceposition));
[ratio_2alpha_alpha_B0Zero_VaryTE_FA5,~,~,mxy_alpha_B0Zero_VaryTE_FA5,mxy_2alpha_B0Zero_VaryTE_FA5,~,mxy_alpha_B0Zero_t0_FA5,mxy_2alpha_B0Zero_t0_FA5]= Mxy_te_Grefz(B1Bias,deg2rad(nominalAngle),kFactor,df_B0Zero,T1,T2,TE,sliceposition,FatSupress);
% FA 65 degrees
nominalAngle = 65;
[ratio_2alpha_alpha_B0Zero_VaryTE_FA65,~,~,mxy_alpha_B0Zero_VaryTE_FA65,mxy_2alpha_B0Zero_VaryTE_FA65,~,mxy_alpha_B0Zero_t0_FA65,mxy_2alpha_B0Zero_t0_FA65]= Mxy_te_Grefz(B1Bias,deg2rad(nominalAngle),kFactor,df_B0Zero,T1,T2,TE,sliceposition,FatSupress);

%======2. Constant Off-resonance: B0=50Hz========
OffResSliceCentre = 50;%v(2);%Hz
nominalAngle = 65;
df_constantB0 = OffResSliceCentre*ones(1,length(sliceposition));
[ratio_2alpha_alpha_constantB0_VaryTE,~,~,mxy_alpha_constantB0_VaryTE,mxy_2alpha_constantB0_VaryTE,~,mxy_alpha_constantB0_t0,mxy_2alpha_constantB0_t0]= Mxy_te_Grefz(B1Bias,deg2rad(nominalAngle),kFactor,df_constantB0,T1,T2,TE,sliceposition,FatSupress);

%% Mxy Magnitude and Phase plots after the RF Excitation - Figure 1 (a), (b)
clc
addpath('DrosteEffect-BrewerMap-ca40391')
% close all
ColorLinePlot = brewermap(9,'Set1');
TE_CentreKspace = 11;
%indexB0Gradient = find(LinearB0gradient_OneSlice==-45);

%----Absolute Mxy Magnitude plots-------
fig=figure();
fig.Color = [1 1 1];
%subplot(2,3,[1 3])
hold on
%t=0
plot(sliceposition,abs(diag(mxy_alpha_B0Zero_t0_FA65)),'LineWidth',2,'Color',ColorLinePlot(3,:));
plot(sliceposition,abs(diag(mxy_2alpha_B0Zero_t0_FA65)),'--','LineWidth',2,'Color',ColorLinePlot(3,:));
plot(sliceposition,abs(diag(mxy_alpha_constantB0_t0)),'LineWidth',2,'Color',ColorLinePlot(4,:));
plot(sliceposition,abs(diag(mxy_2alpha_constantB0_t0)),'--','LineWidth',2,'Color',ColorLinePlot(4,:));

hold off
ylim([0 1.2])
xlim([-1.3 1.3])
xticks(-1.2:0.2:1.2)
xlabel('Z Slice Position (cm)');
ylabel('|Mxy| (a.u.)');
% title('B1+ Factor = 1')
grid on; 
% %Add a), b) or c) to figure
xt = -0.15;
yt = 1.04;
str = {'(a)'};
text(xt,yt,str,'Units','normalized','fontsize',18,'FontName','Arial','FontWeight','bold')
set(gca,'FontSize',18,'FontName','Arial')
legend(    'Mxy_{\alpha} 0 Hz','Mxy_{2\alpha} 0 Hz',...
['M_xy_{\alpha}',num2str(df_constantB0(1,1)),' Hz'],['Mxy_{2\alpha} ',num2str(df_constantB0(1,1)),' Hz'],'Location','NorthEast','FontSize',15,'FontName','Arial','NumColumns',3)
%['Mxy_{\alpha} ',num2str(LinearB0gradient_OneSlice(1,indexB0Gradient)),' Hz/cm'],['Mxy_{2\alpha} ',num2str(LinearB0gradient_OneSlice(1,indexB0Gradient)),' Hz/cm'],...

set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');
figName = 'B0Z_Magnitude';

%----Mxy Phase plots---------

%indexB0Gradient = find(LinearB0gradient_OneSlice==-45);
TE_CentreKspace = 11;%ms
fig=figure();
fig.Color = [1 1 1];
hold on
%t=0
plot(sliceposition,angle(diag(mxy_alpha_B0Zero_t0_FA65)),'LineWidth',2,'Color',ColorLinePlot(3,:));
plot(sliceposition,angle(diag(mxy_2alpha_B0Zero_t0_FA65)),'--','LineWidth',2,'Color',ColorLinePlot(3,:));
plot(sliceposition,angle(diag(mxy_alpha_constantB0_t0)),'LineWidth',2,'Color',ColorLinePlot(4,:));
plot(sliceposition,angle(diag(mxy_2alpha_constantB0_t0)),'--','LineWidth',2,'Color',ColorLinePlot(4,:));
hold off
xlim([-1.3 1.3])
xticks(-1.2:0.2:1.2)
ylim([-pi/2 pi/2])
xlabel('Z Slice Position (cm)');
ylabel('\phi_{Mxy} (radians)');
grid on; 
% title('Mxy phase at time=TE=11ms for B0=0Hz(red), B0=50Hz(blue), dB0/dz=-45Hz/cm(green)')
% %Add a), b) or c) to figure
xt = -0.15;
yt = 1.04;
str = {'(b)'};
text(xt,yt,str,'Units','normalized','fontsize',18,'FontName','Arial','FontWeight','bold')
set(gca,'FontSize',18,'FontName','Arial')
legend(    'Mxy_{\alpha} 0 Hz','Mxy_{2\alpha} 0 Hz',...
['Mxy_{\alpha} ',num2str(df_constantB0(1,1)),' Hz'],['Mxy_{2\alpha} ',num2str(df_constantB0(1,1)),' Hz'],'Location','SouthEast','FontSize',15,'FontName','Arial','NumColumns',3)

set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
figName = 'B0Z_Phase';


%% Mxy Magnitude and Phase plots at time=TE for B0=0Hz, B0=50Hz, Linear B0-gradientZ ±45Hz/cm

% ===== Plot the Results ======
addpath('DrosteEffect-BrewerMap-ca40391')
ColorLinePlot = brewermap(5,'Set1');
TE_CentreKspace = 11;

%Constant Parameters 
T1 = 900; %ms
T2 = 30; %ms
TE = 1:15; %ms %time from the centre of the RF pulse
kFactor = 2; %k*alpha=2alpha
B1Bias = 1;
sliceposition = -1:0.01:1; %cm.
%Negative B0 gradient -45Hz/cm: ratio smaller than on-resonance
B0Gradient=-45;
nominalAngle = 65;
df_LinearB0Gradient = sliceposition.*B0Gradient;
[ratio_2alpha_alpha_NegativeLinearB0Gradient_VaryTE,~,~,mxy_alpha_NegativeLinearB0Gradient_VaryTE,mxy_2alpha_NegativeLinearB0Gradient_VaryTE]= Mxy_te_Grefz(B1Bias,deg2rad(nominalAngle),kFactor,df_LinearB0Gradient,T1,T2,TE,sliceposition,FatSupress);
%Fit the exponential decay to a decaying exponential
fitType = 'a*exp(b*x)';
myFit = fit(TE',ratio_2alpha_alpha_NegativeLinearB0Gradient_VaryTE,fitType,'StartPoint',[1;0.02]);
RatioFactorDecreaseB0Gradient=ratio_2alpha_alpha_NegativeLinearB0Gradient_VaryTE(11,1)/ratio_2alpha_alpha_constantB0_VaryTE(1,1);

fig=figure();
fig.Color = [1 1 1];
hold on
plot(TE,ratio_2alpha_alpha_B0Zero_VaryTE_FA65,'-','LineWidth',2.5,'Color',ColorLinePlot(1,:)); 
plot(TE,ratio_2alpha_alpha_constantB0_VaryTE,'--','LineWidth',2.5,'Color',ColorLinePlot(2,:)); 
plot(TE,ratio_2alpha_alpha_NegativeLinearB0Gradient_VaryTE,'LineWidth',2.5,'Color',ColorLinePlot(3,:)); 
plot([TE_CentreKspace TE_CentreKspace],[ratio_2alpha_alpha_constantB0_VaryTE(1,1) ratio_2alpha_alpha_NegativeLinearB0Gradient_VaryTE(15,1)],'--','LineWidth',2.0,'Color',ColorLinePlot(5,:))
xlabel('TE (ms)');
% ylim([0.6 1.25])
xlim([TE(1) TE(end)])
xticks(1:2:15)
%xticks(TE)
%ylabel('Ratio ($\frac{Mxy_{2\alpha}}{Mxy_{\alpha}}$)','interpreter','latex','FontSize',14,'FontName','Arial');
ylabel('|Mxy_{2\alpha}|/|Mxy_{\alpha}|');
legend('B_0: 0 Hz','Constant B_0: 50 Hz','B_0 Gradient: -45 Hz/cm','Location','SouthWest');
%legend('B0: 0 Hz','Constant B0: 50 Hz','B0-Gradient: -45 Hz/cm','fit to a*exp(-b*x)','Location','SouthWest');
% %Add a), b) or c) to figure
xt = -0.15;
yt = 1.04;
str = {'(a)'};
text(xt,yt,str,'Units','normalized','fontsize',18,'FontName','Arial','FontWeight','bold')

set(gca,'FontSize',18,'FontName','Arial')
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
figName = 'B0Z_RatioTEBelowB0Const';


%Positive B0 gradient 45Hz/cm: ratio larger than on-resonance
B0Gradient=45;
nominalAngle = 65;
df_LinearB0Gradient = sliceposition.*B0Gradient;
[ratio_2alpha_alpha_PositiveLinearB0Gradient_VaryTE,~,~,mxy_alpha_PositiveLinearB0Gradient_VaryTE,mxy_2alpha_PositiveLinearB0Gradient_VaryTE]= Mxy_te_Grefz(B1Bias,deg2rad(nominalAngle),kFactor,df_LinearB0Gradient,T1,T2,TE,sliceposition,FatSupress);
RatioFactorIncreaseB0Gradient=ratio_2alpha_alpha_PositiveLinearB0Gradient_VaryTE(11,1)/ratio_2alpha_alpha_constantB0_VaryTE(1,1);

fig=figure();
fig.Color = [1 1 1];
hold on
plot(TE,ratio_2alpha_alpha_B0Zero_VaryTE_FA65,'-','LineWidth',2.5,'Color',ColorLinePlot(1,:)); 
plot(TE,ratio_2alpha_alpha_constantB0_VaryTE,'--','LineWidth',2.5,'Color',ColorLinePlot(2,:)); 
plot(TE,ratio_2alpha_alpha_PositiveLinearB0Gradient_VaryTE,'LineWidth',2.5,'Color',ColorLinePlot(3,:)); 
% plot(TE_CentreKspace,[ratio_2alpha_alpha_LinearB0Gradient_te(1,1) ratio_2alpha_alpha_B0Zero_VaryTE(11,1)],'Color',ColorLinePlot(5,:),'LineStyle','-')
%plot([TE_CentreKspace TE_CentreKspace],[ratio_2alpha_alpha_constantB0_VaryTE(1,1) ratio_2alpha_alpha_NonLinearB0Gradient_VaryTE(15,1)],'--','LineWidth',2.0,'Color',ColorLinePlot(5,:))
plot([TE_CentreKspace TE_CentreKspace],[0.6 ratio_2alpha_alpha_PositiveLinearB0Gradient_VaryTE(TE_CentreKspace,1)],'--','LineWidth',2.0,'Color',ColorLinePlot(5,:))
hold off
xlabel('TE (ms)');
ylabel('|Mxy_{2\alpha}|/|Mxy_{\alpha}|');
% ylim([0.6 1.25])
xlim([TE(1) TE(end)])
xticks(1:2:15)
% ylabel('Ratio ($\frac{Mxy_{2\alpha}}{Mxy_{\alpha}}$)','interpreter','latex','FontSize',14,'FontName','Arial');
legend('B_0: 0 Hz','Constant B_0: 50 Hz','B_0 Gradient: 45 Hz/cm','Location','SouthWest');
%legend('B0: 0 Hz',['Constant B0: ', num2str(v(2),3),' Hz'],['Positive Non-Linear B0-Gradient: +',num2str(v(1)-v(2),2),' Hz'],'Location','NorthWest');%['B0-Gradient=',num2str(LinearB0gradient_OneSlice(1,indexB0Gradient)),'Hz/cm']);
% %Add a), b) or c) to figure
xt = -0.15;
yt = 1.04;
str = {'(b)'};
text(xt,yt,str,'Units','normalized','fontsize',18,'FontName','Arial','FontWeight','bold')
set(gca,'FontSize',18,'FontName','Arial')
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
figName = 'B0Z_RatioTEAboveB0Const';

%% Figure 2 (a)/(b)

%Explains why the ratio at TE decreases or increases relative to on-resonance
% depending on the B0gradient sign and the slice-select gradient sign. 

% close all
addpath('DrosteEffect-BrewerMap-ca40391')
ColorLinePlot = brewermap(5,'Set1');
fig=figure();
fig.Color = [1 1 1];
tiledlayout(2,1)
nexttile
hold on
plot(sliceposition,abs(mxy_alpha_B0Zero_VaryTE_FA65(:,11)),'LineWidth',2,'Color',ColorLinePlot(3,:));
plot(sliceposition,abs(mxy_2alpha_B0Zero_VaryTE_FA65(:,11)),'--','LineWidth',2,'Color',ColorLinePlot(3,:));
plot(sliceposition,abs(mxy_alpha_NegativeLinearB0Gradient_VaryTE(:,11)),'LineWidth',1.5,'Color',ColorLinePlot(5,:));
plot(sliceposition,abs(mxy_2alpha_NegativeLinearB0Gradient_VaryTE(:,11)),'--','LineWidth',1.5,'Color',ColorLinePlot(5,:));
hold off
xlim([-1 1])
xticks(-1:0.2:1)
xtickangle(0)
ylim([-0.1 1])
xlabel('Z Slice Position (cm)');
ylabel('|Mxy| (a.u.)');
grid on; 
% %Add a), b) or c) to figure
xt = -0.15;
yt = 1.04;
str = {'(a)'};
text(xt,yt,str,'Units','normalized','fontsize',16,'FontName','Arial','FontWeight','bold')
set(gca,'FontSize',16,'FontName','Arial')
nexttile
hold on
plot(sliceposition,angle(mxy_alpha_B0Zero_VaryTE_FA65(:,11)),'LineWidth',1.5,'Color',ColorLinePlot(3,:));
plot(sliceposition,angle(mxy_2alpha_B0Zero_VaryTE_FA65(:,11)),'--','LineWidth',1.5,'Color',ColorLinePlot(3,:));
plot(sliceposition,angle(mxy_alpha_NegativeLinearB0Gradient_VaryTE(:,11)),'LineWidth',1.5,'Color',ColorLinePlot(5,:));
plot(sliceposition,angle(mxy_2alpha_NegativeLinearB0Gradient_VaryTE(:,11)),'--','LineWidth',1.5,'Color',ColorLinePlot(5,:));
plot(sliceposition,angle(mxy_alpha_PositiveLinearB0Gradient_VaryTE(:,11)),'LineWidth',1.5,'Color',ColorLinePlot(2,:));
plot(sliceposition,angle(mxy_2alpha_PositiveLinearB0Gradient_VaryTE(:,11)),'--','LineWidth',1.5,'Color',ColorLinePlot(2,:));
hold off
xlim([-1 1])
xticks(-1:0.2:1)
xtickangle(0)
ylim([-pi pi])
xlabel('Z Slice Position (cm)');
ylabel('\phi_{Mxy} (radians)');
grid on; 
lg=legend('Mxy_\alpha 0 Hz/cm','Mxy_{2\alpha} 0 Hz/cm','Mxy_\alpha -45 Hz/cm','Mxy_{2\alpha} -45 Hz/cm','Mxy_\alpha 45 Hz/cm','Mxy_{2\alpha} 45 Hz/cm','FontSize',14,'FontName','Arial');
lg.Location = 'northeastoutside';
% %Add a), b) or c) to figure
xt = -0.15;
yt = 1.04;
str = {'(b)'};
text(xt,yt,str,'Units','normalized','fontsize',16,'FontName','Arial','FontWeight','bold')

set(gca,'FontSize',16,'FontName','Arial')
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
figName = 'B0Z_WhyRatioDeIncreases';

%% Calculate the B1+ error due to a linear B0 gradient across the slice varying between -45:5:45 Hz/cm

clear;clc;close all

%====== Linear B0-gradientZ ========
LinearB0gradient_OneSlice = -45:5:45;%Hz/cm
nLinearB0Gradient = length(LinearB0gradient_OneSlice);
sliceposition = -1:0.01:1; %cm
nZPositions = length(sliceposition);
nominalAngle = 65;
kFactor = 2; %k*alpha=2alpha
%position between [-1 0 1]cm, with -45Hz/cm you get B0=[-45 0 45] Hz
B1BiasRange = [0.59,1,1.14];
nB1Bias = length(B1BiasRange);

T1 = 900; %ms
T2 = 30; %ms
TE = 1:15; %ms %time from the centre of the RF pulse
FatSupress = 'FS';

%Initialize variables
ratio_2alpha_alpha_VaryLinearB0Gradient_te = zeros(nB1Bias,nLinearB0Gradient);
B1Calculated = zeros(nB1Bias,nLinearB0Gradient);
mxy_alpha_VaryLinearB0Gradient_te = zeros(nZPositions,nLinearB0Gradient,nB1Bias);
mxy_alpha_VaryLinearB0Gradient_t0 = zeros(nZPositions,nLinearB0Gradient,nB1Bias);
mxy_2alpha_VaryLinearB0Gradient_t0 = zeros(nZPositions,nLinearB0Gradient,nB1Bias);
mxy_2alpha_VaryLinearB0Gradient_te = zeros(nZPositions,nLinearB0Gradient,nB1Bias);
ErrorB1 = zeros(nB1Bias,nLinearB0Gradient);
TE_CentreKspace = 11;%ms
OffRes_SliceCentre = 0; %Hz 

%Look-up Surface (LUS) used to calculate the B1+ without B0GradientZ correction and with refocusing gradient 
load('2FASliceProfOffRes100Hz_RatioFSkFactor2_Angles[1_360]deg_Gzref_Position-1_0.1_1.mat','ratioIntenSliceProfOffres','AmbiguityAngleSS','offResonance')
%load('B1MapOptimAngles/2FASliceProfOffRes300Hz_RatioWEkFactor2_Angles[1_360]deg_Position-1_0.1_1.mat','ratioIntenSliceProfOffres','AmbiguityAngleWE','offResonance')
% indexOffres = find(offResonance == round(abs(OffRes_SliceCentre)));%get index of 0Hz off-resonance to interpolate the LUS
% angleAmbiguity = AmbiguityAngleWE(indexOffres,1);%get the ambiguity angle for the 0Hz off-resonance

%----Bloch simulate mxy_alpha and mxy2alpha to calculate the ratio mxy2alpha/mxyalpha with Linear B0-gradientZ-----
for iLinearB0gradient = 1:nLinearB0Gradient    
    for iB1Bias = 1:nB1Bias
        df_LinearB0Gradient = OffRes_SliceCentre+sliceposition.*LinearB0gradient_OneSlice(1,iLinearB0gradient);
        [ratio_2alpha_alpha_LinearB0Gradient_te,integral_mxy_alpha_te,integral_mxy_2alpha_te,mxy_alpha_LinearB0Gradient_te,mxy_2alpha_LinearB0Gradient_te,~,mxy_alpha_LinearB0Gradient_t0,mxy_2alpha_LinearB0Gradient_t0]= Mxy_te_Grefz(B1BiasRange(1,iB1Bias),deg2rad(nominalAngle),kFactor,df_LinearB0Gradient,T1,T2,TE,sliceposition,FatSupress);
        ratio_2alpha_alpha_VaryLinearB0Gradient_te(iB1Bias,iLinearB0gradient) = ratio_2alpha_alpha_LinearB0Gradient_te(TE_CentreKspace,1);
        mxy_alpha_VaryLinearB0Gradient_te(:,iLinearB0gradient,iB1Bias) = mxy_alpha_LinearB0Gradient_te(:,TE_CentreKspace);
        mxy_2alpha_VaryLinearB0Gradient_te(:,iLinearB0gradient,iB1Bias) = mxy_2alpha_LinearB0Gradient_te(:,TE_CentreKspace);
        
        mxy_alpha_VaryLinearB0Gradient_t0(:,iLinearB0gradient,iB1Bias) = diag(mxy_alpha_LinearB0Gradient_t0);
        mxy_2alpha_VaryLinearB0Gradient_t0(:,iLinearB0gradient,iB1Bias) = diag(mxy_2alpha_LinearB0Gradient_t0);
        
        
        % Estimate B1+ WITH GradientB0Z Correction
        %B1_est = lsqnonlin(@(B1Bias)  Mxy_te_gradientB0SliceProfile(B1Bias,deg2rad(nominalAngle),kFactor,OffRes_SliceCentre,T1,T2,TE_CentreKspace,sliceposition,SliceThicknessFactor,FatSupress)-ratio_2alpha_alpha_VaryLinearB0Gradient_te(iB1Bias,iLinearB0gradient),B1Guess,[],[],opts);

        % Estimate B1+ WITHOUT GradientB0Z Correction
        %FS
        angleEstimated = interp1(ratioIntenSliceProfOffres(1:AmbiguityAngleSS(1),101),1:AmbiguityAngleSS(1),ratio_2alpha_alpha_VaryLinearB0Gradient_te(iB1Bias,iLinearB0gradient),'linear',NaN);
        %WE
        %angleEstimated = interp1(ratioIntenSliceProfOffres(1:angleAmbiguity,indexOffres),1:angleAmbiguity,ratio_2alpha_alpha_VaryLinearB0Gradient_te(iB1Bias,iLinearB0gradient),'linear',NaN);
        B1Calculated(iB1Bias,iLinearB0gradient) = angleEstimated./nominalAngle;
        ErrorB1(iB1Bias,iLinearB0gradient) = (B1Calculated(iB1Bias,iLinearB0gradient) - B1BiasRange(1,iB1Bias));%/B1BiasRange(1,iB1Bias);

    end
end

% ======== Figure 2 (c) ===============
%Linear fit B1+ error=f(dB0/dz)
indexZeroB0Gradient = find(LinearB0gradient_OneSlice==0);
indexB1Factor1 = find(B1BiasRange==1);
ratioZeroB0Gradient_B1Factor1 = ratio_2alpha_alpha_VaryLinearB0Gradient_te(indexB1Factor1,indexZeroB0Gradient);

[p,S] = polyfit(LinearB0gradient_OneSlice,ErrorB1(indexB1Factor1,:),1);
[y_fit,delta] = polyval(p,LinearB0gradient_OneSlice,S);
RSquared = 1 - (S.normr/norm(ErrorB1(indexB1Factor1,:) - mean(ErrorB1(indexB1Factor1,:))))^2;
close all
ColorLinePlotSet1 = brewermap(9,'Set1');
fig=figure();
fig.Color = [1 1 1];
hold on
%B1Factor=0.59
plot(LinearB0gradient_OneSlice,ErrorB1(1,:),'-.*','MarkerSize',9,'LineWidth',2,'Color',ColorLinePlotSet1(5,:)); 
%B1Factor=1.0
plot(LinearB0gradient_OneSlice,ErrorB1(indexB1Factor1,:),'-.*','MarkerSize',9,'LineWidth',2,'Color',ColorLinePlotSet1(3,:)); 
%B1Factor=1.14
plot(LinearB0gradient_OneSlice,ErrorB1(nB1Bias,:),'-.*','MarkerSize',9,'LineWidth',2,'Color',ColorLinePlotSet1(2,:)); 
%plot linear fit y_fit and 95% prediction interval y±2?
plot(LinearB0gradient_OneSlice,y_fit,'--','Color',ColorLinePlotSet1(8,:),'LineWidth',2)
% plot(LinearB0gradient_OneSlice,y_fit+2*delta,'m--',LinearB0gradient_OneSlice,y_fit-2*delta,'m--')
line([LinearB0gradient_OneSlice(1) LinearB0gradient_OneSlice(end)],[0 0],'Color','k','LineStyle','-.')
hold off
xticks(-45:15:45)
ylim([-0.25 0.25])
xlabel('B_0 gradient through slice (Hz/cm)');
ylabel('B_1^+ Factor Error');
legend('B_1^+ Factor = 0.59','B_1^+ Factor = 1.0','B_1^+ Factor = 1.14','Linear Fit to B_1^+ = 1.0','Location','NorthEast')
set(gca,'FontSize',12,'FontName','Arial')
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');

% str = {'Linear Fit: B_1^+ Error=f(dB_0/dz)',...
%     ['y = ',num2str(p(1),2), 'x',num2str(p(2),2)],['R^{2} = ',num2str(RSquared,4)]};
% dim = [0.2 .06 0.11 0.22]; %Specify dim as a four-element vector of the form [x y w h].
% annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontName','Arial');
figName = 'B0Z_B1Error';
