function [ratio_te,integral_mxy_te_alpha,integral_mxy_te_2alpha,mxy_VaryTE,mxy_2alpha_VaryTE,ratio_t0,mxy_t0,mxy_t0_2alpha] = Mxy_te_Grefz(B1,nominalAngle,kFactor,df,T1,T2,TE,sliceposition,FatSup)

%Date: 04122020
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2023

%Summary: calculate the Signal130/Sginal65 assuming varying B0 across the slice profile
%Input: 
    %B1 : unkown parameter to estimate
    %nominalAngle: alpha (65 degrees)
    %kfactor: k*alpha (2*alpha)
    %df: off-resonance frequency vector
    %T1: relaxation time
    %T2: decay
    %TE: echo time
    %sliceposition: z position for Bloch simulation
    %FatSup: Fat Saturation

%Output: 
    %ratio_te: ratio of integrated signals Signal2alpha/Signalalpha for a range of TEs
    %integral_mxy_te_alpha:integrated signal Signalalpha for a range of TEs
    %integral_mxy_te_2alpha:integrated signal Signal2alpha for a range of TEs
    %mxy_VaryTE: mxy complex signal alpha for a range of TEs
    %mxy_2alpha_VaryTE: mxy complex signal 2alpha for a range of TEs
    %ratio_t0: ratio of integrated signals Signal2alpha/Signalalpha immediately after excitation pulse
    %mxy_t0: mxy complex signal alpha immediately after excitation pulse
    %mxy_t0_2alpha: mxy complex signal 2alpha immediately after excitation pulse

%path to Bloch simulation: Wind_blochMxyRFInterp_gradientB0SliceProfile
%%IMP: add path here to Prof. Hargreaves bloch simulator
%path to Bloch simulation: Wind_blochMxyRFInterp_gradientB0SliceProfile

%Note: it does not matter what T1 and T2 you use in the simulation
%(Wind_blochMxyRFInterp_gradientB0SliceProfile) to calculate the signal at end of RF
%pulse as there should be no T1 or T2 relaxation during the time the RF
%pulse is played out. Confirmed this with T1 and T2 of:
% T1 = 1000; %ms
% T2 = 1000; %ms
% gives the same answer as with T1 = 800ms and T2=40ms.
SliceThicknessFactor = 1;
if strcmp(FatSup(1:2),'FS')
    [mxy_t0,mx,my,mz] = Wind_blochMxyRFInterp_gradientB0SliceProfile_Gzref(B1,nominalAngle,df,sliceposition,SliceThicknessFactor,T1,T2);
    [mxy_t0_2alpha,mx_2alpha,my_2alpha,mz_2alpha] = Wind_blochMxyRFInterp_gradientB0SliceProfile_Gzref(B1,nominalAngle*kFactor,df,sliceposition,SliceThicknessFactor,T1,T2);
elseif strcmp(FatSup(1:2),'WE')
    [mxy_t0,mx,my,mz] = WE_blochMxyRFInterpGradShift(B1,nominalAngle,df,sliceposition,T1,T2);
    [mxy_t0_2alpha,mx_2alpha,my_2alpha,mz_2alpha] = WE_blochMxyRFInterpGradShift(B1,nominalAngle*kFactor,df,sliceposition,T1,T2);   
end
%Calculate the Ratio 130/65 with mxy signal at the end of the excitation RF pulse: t=16870 microseconds
[ratio_t0] = ratio_kAlpha_Alpha(mxy_t0,mxy_t0_2alpha);

%Calulate the Mxy for each slice position and T1 relaxation and T2 decay
%after each time TE. It uses the signal calculated at the end of the RF
%pulse. Each slice position has a certain off-resonance frequency saved in
%df which will result in a different propagation matrix
nsamples_sliceProfile = length(sliceposition);
mxy_VaryTE = zeros(nsamples_sliceProfile,length(TE));
mxy_2alpha_VaryTE = zeros(nsamples_sliceProfile,length(TE));
integral_mxy_te_alpha = zeros(length(TE),1);
integral_mxy_te_2alpha = zeros(length(TE),1);


Mx_te = zeros(nsamples_sliceProfile,1);
My_te = zeros(nsamples_sliceProfile,1);
Mz_te = zeros(nsamples_sliceProfile,1);
Mx_te_2alpha = zeros(nsamples_sliceProfile,1);
My_te_2alpha = zeros(nsamples_sliceProfile,1);
Mz_te_2alpha = zeros(nsamples_sliceProfile,1);

%get the time (microseconds) where the peak of the RF pulse occurs
% Interpolate RF Pulse and Gradient Z to 1 microsend time
% interval
%[value,loc] = max(rf_phase_interp(:));
%PeakRFPulse = time_interp (1,loc) %15650, 17470
%EndGrzRef = time_interp (1,end) %microseconds
%EndGrzRef-PeakRFPulse = 1820 %microseconds
if strcmp(FatSup(1:2),'FS')
    Time_CentreRFUntilGzRefEnd = (1820)*10^-3; %Gzref = 22.8mT/m
    %Time_CentreRFUntilGzRefEnd = (2443)*10^-3; %Gz ref = 5.96mT/m
    freePrecessionTime = TE -Time_CentreRFUntilGzRefEnd;% TE is defined from the centre of the RF pulse, but free precession only starts after the Gz Ref ends
elseif strcmp(FatSup(1:2),'WE') 
    Time_CentreRFUntilGzRefEnd = (2887.5)*10^-3; %ms
    freePrecessionTime = TE -Time_CentreRFUntilGzRefEnd;% TE is defined from the centre of the RF pulse, but free precession only starts after the Gz Ref ends
end

for iTE = 1:length(freePrecessionTime) %loop over different echo times
    for islicePos = 1:nsamples_sliceProfile %loop over the slice profile

        M_alpha = [mx(islicePos,islicePos);my(islicePos,islicePos);mz(islicePos,islicePos)]; %Magnetization immediately after you play excitation RF pulse alpha
        M_2alpha = [mx_2alpha(islicePos,islicePos);my_2alpha(islicePos,islicePos);mz_2alpha(islicePos,islicePos)]; %Magnetization immediately after you play excitation RF pulse 2*alpha
        
        % ===== Get the Propagation Matrix ======
        %propagates the signal from end of RF pulse until time TE
        [Ate,Bte] = freeprecess(freePrecessionTime(iTE),T1,T2,df(1,islicePos));
        
        
        Mte = Ate*M_alpha+Bte;	% Magnetization at t=TE for alpha
        Mte_2alpha = Ate*M_2alpha+Bte;	% Magnetization at t=TE for 2*alpha
        
        
        Mx_te(islicePos,1) = Mte(1,1);%Mx is the first component
        My_te(islicePos,1) = Mte(2,1);%My is the second component
        Mz_te(islicePos,1) = Mte(3,1);%Mz is the third component
        
        Mx_te_2alpha(islicePos,1) = Mte_2alpha(1,1);
        My_te_2alpha(islicePos,1) = Mte_2alpha(2,1);
        Mz_te_2alpha(islicePos,1) = Mte_2alpha(3,1);
    end
    
    mxy_te = Mx_te+1i.*My_te; %mxy for alpha
    mxy_VaryTE(:,iTE) = mxy_te;
    mxy_te_2alpha = Mx_te_2alpha+1i.*My_te_2alpha; %mxy for 2*alpha
    mxy_2alpha_VaryTE(:,iTE) = mxy_te_2alpha;
        
    integral_mxy_te_alpha(iTE,1) = abs(trapz(mxy_te(:,1))); %integrate the signal across the slice profile
    integral_mxy_te_2alpha(iTE,1) = abs(trapz(mxy_te_2alpha(:,1)));
    
end

%Ratio of integrated Signal2alpha and integrated Signalalpha
ratio_te = integral_mxy_te_2alpha./integral_mxy_te_alpha;



end

