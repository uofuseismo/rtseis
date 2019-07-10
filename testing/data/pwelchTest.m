function pwelchTest()
    fs = 10000;  % Sampling rate Hz
    freq = 1234; % Frequency of sine wave
    N = 10000;
    time = linspace(0, N-1, N)/fs; % Time
    x = sin(2*pi*freq*time); % Evaluate sine wave at times
    % A traditional test
    w = hann(1024);
    overlap = 512;
    fftLength = length(w);
    [pxxPower, f] = pwelch(x, w, overlap, fftLength, fs, 'power');
    [pxxPSD,   f] = pwelch(x, w, overlap, fftLength, fs, 'psd');
    file = fopen('welchTest1.txt', 'w');
    for i=1:length(f)
        fprintf(file, '%.10e %.10e %.10e\n', f(i), pxxPower(i), pxxPSD(i));
    end
    fclose(file);
    % This is a strange test intended to test the non-default settings
    w = hamming(1013);
    overlap = 501;
    fftLength = 1091;
    [pxxPower, f] = pwelch(x, w, overlap, fftLength, fs, 'power');
    [pxxPSD,   f] = pwelch(x, w, overlap, fftLength, fs, 'psd');
    file = fopen('welchTest2.txt', 'w');
    for i=1:length(f)
        fprintf(file, '%.10e %.10e %.10e\n', f(i), pxxPower(i), pxxPSD(i)); 
    end 
    fclose(file);
end 
