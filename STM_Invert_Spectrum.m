function [Recon, newTF,TFreference] = STM_Invert_Spectrum(MSstruct,initialPhase,iter)
%
%   Invert the spectrum from a filtered TF (spectrogram) to a waveform
% [Recon, newTF,TFreference] = STM_Invert_Spectrum(MSstruct,initialPhase,iter)
%
debug = 0;

TF = MSstruct.new_TF;       % Time frequency is the new (filtered) TF outputed by STM_Filter_Mod

if MSstruct.Args.MS_log      % Revert to analytic amplitude if power is in log (DB) units
    mx = max(max(TF));
    TF(find(TF<(0.05*mx))) = 0;
    TF = realpow(10.0,(TF-MSstruct.TF_Args.DBNoise)./20.0);
    TF = TF-min(min(TF));
else
    TF = sqrt(TF);                   % Revert to analytic amplitude (from power)
end

if ~exist('iter'), iter = 20;
end
if  ~exist('initialPhase')
    initialPhase =  [];
end


switch MSstruct.TF_Args.CB_Filter
    case 'gauss'
        nchans = MSstruct.TF_Args.Nchans;
        [p,q] = rat(MSstruct.TF_Args.TF_ReFs /MSstruct.TF_Args.Fs);
        for ch = 1:nchans
                %clear thisBand;
                TFreference(ch,:) = resample(TF(ch,:),q,p);                      % upsample
        end
        if isempty(initialPhase) 
            initialPhase = (rand(size(TFreference))-0.5*ones(size(TFreference))).*2*pi;
        end
        TFreference = TFreference(:,1:size(initialPhase,2));
        %TFreference = TFreference./max(max(TFreference));
        newTF = TFreference.*exp(j.*initialPhase);
        
        for i = 1:iter    
            % Invert spectrogram
            thisSignal = zeros(1,size(newTF,2));
            for ch=1:nchans
                clear thisTF thisBand;
                %[p1,q1] = rat(MSstruct.TF_Args.Fs /round(MSstruct.TF_Args.CB_FrHigh(ch)*2));
                %thisTF = resample(newTF(ch,:),p1,q1); %downsample
                %thisBand = resample(inverseHilbert(thisTF,smooth_F),q1,p1); % inverse and upsample
                %thisBand(isnan(thisBand)) = 0;
                %thisBand = thisBand(1:length(thisSignal));
                %if sum(isnan(thisBand)), keyboard; end
               % if ch == 40, keyboard; end
                thisBand = inverseHilbert(newTF(ch,:));
                
                thisSignal = thisSignal + band_pass(thisBand,MSstruct.TF_Args.Fs,MSstruct.TF_Args.CB_FrLow(ch),MSstruct.TF_Args.CB_FrHigh(ch),1,'HMFWgauss');                                       % sum bands into one signal
               %thisSignal = thisSignal + thisBand;                                       % sum bands into one signal
            end
            % Calculate new spectrogram
            clear newTF*;
            for ch = 1:nchans
                newTF(ch,:) = hilbert(band_pass(thisSignal,MSstruct.TF_Args.Fs,MSstruct.TF_Args.CB_FrLow(ch),MSstruct.TF_Args.CB_FrHigh(ch),1,'HMFWgauss'));
                %newTFdn(ch,:) = resample(abs(newTF(ch,:)),p,q); 
            end
            %newTFdn = newTFdn +abs(min(min(newTFdn)));
            %newTFsv = abs(newTF);
            %newTF = newTF./max(newTFsv(:)).*max(TFreference(:));
            theErr = sum(sum( ...
                               (abs(newTF)./max(max(abs(newTF)))-TFreference./max(max(TFreference))).^2))/ ...
                          sum(sum( (TFreference./max(max(TFreference))).^2 ))*100;
            
            %theErr = sum(( abs(newTF(:)) -abs(TFreference(:)) ).^2) ./ sum(abs(TFreference(:)))*100;
            
            if debug,
                subplot(311);
                imagesc(abs(TFreference)); axis xy;colorbar;
                subplot(312);
                imagesc(abs(newTF)); colorbar;axis xy; title('Recon');
                subplot(313); plot(thisSignal);colorbar; xlim(size(thisSignal)); drawnow;
                %sound(thisSignal,MSstruct.TF_Args.Fs);
                set(gcf,'color','w');
                % quantify and output difference
                fprintf('Iteration %d; %2.3f difference between reference and new spectrogam\n',i,theErr);

%                keyboard;
            else
                fprintf('Iteration %d; %2.3f difference between reference and new spectrogam\n',i,theErr);
            end

            newTF = TFreference.*exp(j.*angle(newTF)); % use magnitude from the filtered TF and phase from the new TF

            
        end
          Recon = thisSignal;
    case 'FIR'
        error('FIR Spectrum inversion not yet implemented');
    case 'Shamma'
        try
            loadload
         catch ME
             error('nsltools functions not found');
        end
        Recon = aud2wav(TF',[],[MSstruct.TF_Args.CB_paras iter 0 0]);
end
            
end
function out = inverseHilbert(band)
    
    out = zeros(size(band));
     thisband = band;
        n = length(thisband);
        B = fft(thisband);
        
        if ~mod(n,2) % even
            B([2:n/2 n/2+2:end]) = B([2:n/2 n/2+2:end])./2;
            %B(n/2+2:end) = conj(B(n/2:-1:2));
        else                % odd
            B(2:end) = B(2:end)./2;
            %B(ceil(n/2)+1:end) = conj(B(ceil(n/2):-1:2));
        end
     out = ifft(B,'symmetric');
end

        