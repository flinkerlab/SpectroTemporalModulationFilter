function [MSout] = STM_Filter_Mod(TFstruct,TemModFilter, SpecModFilter,Args,Type)
%
%  Caluclate modulation spectrum and filter
% 
%   MSout = STM_Filter_Mod(TFstruct, <TemModFilter, SpecModFilter,Args,Type>)
%           TFstruct                Time frequency representation from STM_CreateTF
%           TemModFilter            Temporal Modulation filter cutoff
%           SpecModFilter           Spectral Modulation filter cutoff
%           Args                    Arguments, defaults given by Args = STM_Filter_Mod
%           Type                    default 'lowpass'. supported 'highoass', 'bandpass'
%
%           MSout                       output MS structure
%               .new_TF                New TF after MS filtering
%               .new_MS               New MS amplitude after filtering
%               .new_MS_phase   New MS phase after filtering
%               .orig_TF                Original TF 
%               .orig_MS               Original MS amplitude (no filtering)
%               .orig_MS_phase   Original MS phase (no filtering)
%               .x_axis                   MS x axis labels for plotting
%               .y_axis                   MS y axis labels for plotting
%               .Args                     Arguments used to create the MS structure
%               .Args_TF               Arguments used to create the TF  structure (passed from STM_CreateTF)
%
%   Examples:
%       TF = STM_CreateTF(signal,16000); % calculate the TF matrix
%       MS = STM_Filter_Mod(TF);
%       imagesc(MS.x_axis,MS.y_axis,log(MS.orig_MS)); axis xy;
%       imagesc(MS.x_axis,MS.y_axis(fix(length(MS.y_axis)/2+1):end),log(MS.orig_MS(fix(length(MS.y_axis)/2+1):end,:))); axis xy;
%     
%       MS = STM_Filter_Mod(TF,10,15);
%       imagesc(MS.new_TF); axis xy;
%       
%
%       Adeen Flinker, Jan 2013 (adeen.f@gmail.com)
%

if ~exist('Args') || isempty(Args)
    Args.MS_log = 1;                 % work on log transformed spectrogram in DB (1) or analytic amplitude (0)
    Args.MS_filter = 'lowpass'; % lowpass, highpass or bandpass 
    Args.MS_fstep       = 1/28;        % Default 29 chans per ocatve (default in TFstruct.Args.CF_CenterFrs)
    Args.MS_keep_positive_negative_both = 'both'; %default is to keep both positive and negative modulations. positive/negative is supported
end

if exist('Type')
    Args.MS_filter = Type;
end

if strcmp(Args.MS_filter,'highpass')
    filt_dir = -1;
else
    filt_dir = 1;
end
if ~exist('SpecModFilter') 
    SpecModFilter = filt_dir*Inf;
elseif isempty(SpecModFilter)
    SpecModFilter = filt_dir*Inf;
end
if ~exist('TemModFilter') 
    TemModFilter = filt_dir*Inf;
elseif isempty(TemModFilter)
    TemModFilter = filt_dir*Inf;
end


if nargin<1, MSout = Args; return; end

if Args.MS_log
    TF = TFstruct.TFlog;
else
    TF = TFstruct.TF;
end

orig_ms = fft2(TF);

[Nfr,Ntm] = size(TF);
fstep = Args.MS_fstep;
ampsrate = TFstruct.Args.TF_ReFs;

for i=1:ceil((Nfr+1)/2)
    dwf(i)= (i-1)*(1/(fstep*Nfr));    % positive spectral modulation frequencies
    if (i > 1)
        dwf(Nfr-i+2)=-dwf(i);          % negative spectral modulation frequencies
    end
end
for i=1:ceil((Ntm+1)/2)
    dwt(i) = (i-1)*(ampsrate/Ntm);
    if (i > 1 )
        dwt(Ntm-i+2) = -dwt(i);
    end
end

dfi=0.16;%0.0001;                           % 0.0001 cycle per octave ramp in frequency
dti=1;                                               % One Hz ramp in time

newamp = abs(orig_ms);
newphase = angle(orig_ms);
gainmap=ones(size(orig_ms)); % Define a gain by which to multiply the mod spectrum

if TemModFilter == Inf & SpecModFilter == Inf
        Args.MS_filter = 'do_nothing';
end

% code below modified from Eliott et al. modfilter.m
% Sepcifically, the cosine ramp and gain definition for filtering the
% modulation domain ( for f=1:Nfr, for t=1:Ntm, cosine ramp)
% notch, bandpass and lowpass were modified to include proper ramp
% definition and take into accound different algorithmic end cases.
%
%
% This defines the ramp of the gain from 0 to 1
switch  Args.MS_filter
   case 'notch'
       if length(SpecModFilter)<2 || isempty(SpecModFilter)
           wf_high = -Inf;
           wf_low = Inf;
       else
        wf_high = SpecModFilter(2);
        wf_low = SpecModFilter(1);           
       end
       if length(TemModFilter)<2 || isempty(TemModFilter)
           wt_high = -Inf;
           wt_low = Inf;
       else
        wt_high = TemModFilter(2);
        wt_low = TemModFilter(1);   
       end
        
        for f=1:Nfr
         for t=1:Ntm
            % Define the regions to set to zero gain - first along the
            % wf axis
            if ((abs(dwf(f)))>wf_low) && ((abs(dwf(f)))<wf_high)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
                newphase(1,1) = 0.0;  % The phase of the DC value has to be zero
            end

            if (wf_high~=0)
                if ((abs(dwf(f))>=wf_low-dfi) && (abs(dwf(f))<=(wf_low))) || ...
                   ((abs(dwf(f))>=wf_high) && (abs(dwf(f))<=(wf_high+dfi)))
                    gainmap(f,t)=gainmap(f,t)*sin((((abs(dwf(f))-wf_high)./dfi))*(pi/2))^2;
                end
            end

            % Define the regions to set to zero - along the wt axis
            if ((abs(dwt(t)))>wt_low) && ((abs(dwt(t))) <wt_high)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
            end

            if (wt_high~=0)
                if ((abs(dwt(t))>=(wt_low-dti) && (abs(dwt(t))<=(wt_low)))) || ...
                    ((abs(dwt(t))>=(wt_high) && (abs(dwt(t))<=(wt_high+dti))))
                    gainmap(f,t)=gainmap(f,t)*sin((((abs(dwt(t)))-wt_high./dti))*(pi/2))^2;
                end
            end
        end
    end

        case 'bandpass'
       if length(SpecModFilter)<2 || isempty(SpecModFilter)
           wf_high = Inf;
           wf_low = -Inf;
       else
        wf_high = SpecModFilter(2);
        wf_low = SpecModFilter(1);           
       end
       if length(TemModFilter)<2 || isempty(TemModFilter)
           wt_high = Inf;
           wt_low = -Inf;
       else
        wt_high = TemModFilter(2);
        wt_low = TemModFilter(1);   
       end

      
        for f=1:Nfr
         for t=1:Ntm
            % Define the regions to set to zero gain - first along the
            % wf axis
            if ((abs(dwf(f)))>wf_high+dfi) || ((abs(dwf(f)))<wf_low-dfi)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
                newphase(1,1) = 0.0;  % The phase of the DC value has to be zero
            end

            if (wf_high~=0)
                if ((abs(dwf(f))>=wf_low-dfi) && (abs(dwf(f))<=(wf_low))) || ...
                   ((abs(dwf(f))>=wf_high) && (abs(dwf(f))<=(wf_high+dfi)))
                    gainmap(f,t)=gainmap(f,t)*cos((((abs(dwf(f))-wf_high)./dfi))*(pi/2))^2;
                end
            end

            % Define the regions to set to zero - along the wt axis
            if ((abs(dwt(t)))>wt_high+dti) || ((abs(dwt(t))) <wt_low-dti)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
            end

            if (wt_high~=0)
                if ((abs(dwt(t))>=(wt_low-dti) && (abs(dwt(t))<=(wt_low)))) || ...
                    ((abs(dwt(t))>=(wt_high) && (abs(dwt(t))<=(wt_high+dti))))
                        
                    gainmap(f,t)=gainmap(f,t)*cos((((abs(dwt(t))-wt_high)./dti))*(pi/2))^2;
                end
            end
        end
    end

    case 'highpass'
        
        if ~exist('SpecModFilter') 
            SpecModFilter = -Inf;
        elseif isempty(SpecModFilter)
            SpecModFilter = -Inf;
        end
        if ~exist('TemModFilter') 
            TemModFilter = -Inf;
        elseif isempty(TemModFilter)
            TemModFilter = -Inf;
        end

        wf_high = SpecModFilter;
        wt_high = TemModFilter;
        
        for f=1:Nfr
         for t=1:Ntm
            % Define the regions to set to zero gain - first along the
            % wf axis
            if ((abs(dwf(f)))<wf_high)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
               % newphase(1,1) = 0.0;  % The phase of the DC value has to be zero
            end

            if (wf_high~=0)
                if ((abs(dwf(f))>=wf_high) && (abs(dwf(f))<=(wf_high+dfi)))
                    gainmap(f,t)=gainmap(f,t)*cos((((abs(dwf(f))-(wf_high-dfi))./dfi))*(pi/2))^2;
                end
            end

            % Define the regions to set to zero - along the wt axis
            if ((abs(dwt(t)))<wt_high)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
            end

            if (wt_high~=0)
                if ((abs(dwt(t))>=(wt_high) && (abs(dwt(t))<=(wt_high+dti))))
                    gainmap(f,t)=gainmap(f,t)*cos((((abs(dwt(t))-(wt_high-dti))./dti))*(pi/2))^2;
                end
            end
        end
    end
    
    case 'lowpass'
         if ~exist('SpecModFilter') 
            SpecModFilter = Inf;
        elseif isempty(SpecModFilter)
            SpecModFilter = Inf;
        end
        if ~exist('TemModFilter') 
            TemModFilter = Inf;
        elseif isempty(TemModFilter)
            TemModFilter = Inf;
        end
        wf_high = SpecModFilter;
        wt_high = TemModFilter;
        
        for f=1:Nfr
         for t=1:Ntm
            % Define the regions to set to zero gain - first along the
            % wf axis
            if ((abs(dwf(f)))>wf_high+dfi)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
                newphase(1,1) = 0.0;  % The phase of the DC value has to be zero
            end

            if (wf_high~=0)
                if ((abs(dwf(f))>=wf_high) && (abs(dwf(f))<=(wf_high+dfi)))
                    gainmap(f,t)=gainmap(f,t)*cos((((abs(dwf(f))-wf_high)./dfi))*(pi/2))^2;
                end
            end

            % Define the regions to set to zero - along the wt axis
            if ((abs(dwt(t)))>wt_high+dti)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
            end

            if (wt_high~=0)
                if ((abs(dwt(t))>=(wt_high) && (abs(dwt(t))<=(wt_high+dti))))
                    gainmap(f,t)=gainmap(f,t)*cos((((abs(dwt(t))-wt_high)./dti))*(pi/2))^2;
                end
            end
        end
    end
    
end
    newphase(1,1) = 0.0;  % The phase of the DC value has to be zero

if isfield(Args,'MS_keep_positive_negative_both')
   switch Args.MS_keep_positive_negative_both
       case 'both',                              % do nothing
%        case 'positive',gainmap(:,dwt<0) = 0; 
%                         newphase(:,dwt<0)=(rand(size(newphase,1),sum(dwt<0))-0.5)*2*pi;
%                        
%        case 'negative',gainmap(:,dwt>0) = 0; 
%                         newphase(:,dwt>0)=(rand(size(newphase,1),sum(dwt>0))-0.5)*2*pi;
                       
        case 'positive',gainmap(dwf>0,dwt<0) = 0;
                        gainmap(dwf<0,dwt>0) = 0; 
                        newphase(dwf>0,dwt<0)=(rand(sum(dwf>0),sum(dwt<0))-0.5)*2*pi;
                        newphase(dwf<0,dwt>0)=(rand(sum(dwf<0),sum(dwt>0))-0.5)*2*pi;
        case 'negative',gainmap(dwf>0,dwt>0) = 0;
                        gainmap(dwf<0,dwt<0) = 0; 
                        newphase(dwf>0,dwt>0)=(rand(sum(dwf>0),sum(dwt>0))-0.5)*2*pi;
                        newphase(dwf<0,dwt<0)=(rand(sum(dwf<0),sum(dwt<0))-0.5)*2*pi;
                         %newphase(:,dwt>0)=(rand(size(newphase,1),sum(dwt>0))-0.5)*2*pi;
                         %newphase(dwf>0,:)=(rand(sum(dwf>0),size(newphase,2))-0.5)*2*pi;
       otherwise, error('unkown MS_keep_positive_negative_both value. Please use ''both'' ''positive'' ''negative''.');
end
newamp = newamp.*gainmap;
new_ms = newamp.*exp(complex(0,newphase));

MSout.new_TF = real(ifft2(new_ms));
MSout.new_MS = fftshift(abs(new_ms));
MSout.new_MS_phase =fftshift(angle(new_ms));
MSout.orig_TF = TF;
MSout.orig_MS = fftshift(abs(orig_ms));
MSout.orig_MS_phase = fftshift(angle(orig_ms));
MSout.x_axis = fftshift(dwt);
MSout.y_axis = fftshift(dwf);
MSout.Args = Args;
MSout.TF_Args = TFstruct.Args;

end

