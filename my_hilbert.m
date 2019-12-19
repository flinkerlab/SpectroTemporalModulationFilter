function [filt_signal]=my_hilbert(input, sampling_rate, lower_bound, upper_bound, tm_OR_fr,wind)
%  [filt_signal]=my_hilbert(input, sampling_rate, lower_bound, upper_bound, tm_OR_fr)
%     input         - input signal to be filtered (time or frequency domain)
%     sampling_rate - signal's sampling rate
%     lower_bound     - lower frequency bound for bandpass filtering
%     upper_bound   - upper frequency bound for bandpass filtering
%     tm_OR_fr      - 1 if the input signal is in the time domain, 0 if it
%                     is in the frequency domain
%     wind          - windows type:
%                       'HMFWgauss'     - HMFW of upper_bound-lower_bound
%                       'flatgauss' - gaussian with a the maximum point
%                                     flatened to upper_bound-lower_bound
%                                     length
%
%  The function returns the filtered hilbert signal (low->high) in the time domain
%  Written by Adeen Flinker and Lavi Secundo (circa 2008)
%
    if (nargin<6)
        wind = 'flatgauss';
    end
    if (nargin<5)
        tm_OR_fr=1;
    end
    if (nargin<4)
        error('Please enter at least 4 arguments');
    end
    
    if (size(input,1)==1 || size(input,2)==1) && (size(input,2)==1)
       input = input';
    end
    max_freq=sampling_rate/2;
    df=2*max_freq/length(input);
    centre_freq=(upper_bound+lower_bound)/2;
    filter_width=upper_bound-lower_bound;
    x=0:df:max_freq;
      gauss_width = 1;

    if (isnumeric(wind))
        gauss_width = wind;
        wind = 'flatgauss';
    end
    
    switch (wind)
        case 'flatgauss',
                        gauss=exp(-(x-centre_freq).^2.*gauss_width);
                        cnt_gauss = round(centre_freq/df);
                        flat_padd = round(filter_width/df);  % flat padding at the max value of the gaussian
                        padd_left = floor(flat_padd/2);
                        padd_right = ceil(flat_padd/2); 
                    	our_wind = [gauss((padd_left+1):cnt_gauss) ones(1,flat_padd) gauss((cnt_gauss+1):(end-padd_right))];
        case 'HMFWgauss',
                        sigma = filter_width./(2*sqrt(2*log(2)));          % standrad deviation to conform with HMFW of filter_width
                        gauss=exp(-(x-centre_freq).^2./(2*sigma^2));
                        our_wind = gauss;
        otherwise,      error('no valid window');
     end	
    our_wind=[our_wind zeros(1,length(input)-length(our_wind))];
    if (tm_OR_fr==1)
        input=fft(input,[],2);
    end       
    our_wind(1) = our_wind(1)/2; % DC component is halved
    our_wind = repmat(our_wind,size(input,1),1).*2;
    filt_signal=ifft(input.*our_wind,[],2);

end
    

