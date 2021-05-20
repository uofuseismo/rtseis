function spectrogramTest()
   nSamp = 2048;
   Fs = 1024;
   t = (0:nSamp-1)'/Fs;
   t1 = t(1:nSamp/2);

   x11 = sin(2*pi*400*t1);
   x12 = chirp(t1-t1(nSamp/4),150,nSamp/Fs,1750,'quadratic');
   x1 = x11+x12;
   
   t2 = t(nSamp/2+1:nSamp);

   x21 = chirp(t2,400,nSamp/Fs,100);
   x22 = chirp(t2,550,nSamp/Fs,250);
   x2 = x21+x22;
   
   SNR = 20;
   rng('default')

   sig = [x1;x2];
   sig = sig + randn(size(sig))*std(sig)/db2mag(SNR);
   
   fileID = fopen("spectrogramNoisyChirpSignal.txt", "w");
   for i =1:length(sig)
       fprintf(fileID, "%.15e\n", sig(i));
   end
   fclose(fileID);
   
   nwin = 63;
   wind = kaiser(nwin,17);
   nlap = nwin-10;
   nfft = 256;
   
   [S,F,T] = spectrogram(sig, wind, nlap, nfft, Fs);
   %F
   %T
   sShape = size(S);
   fileID = fopen("spectrogramNoisyChirp.txt", "w");
   for i=1:sShape(1)
       for j=1:sShape(2)
           fprintf(fileID, "%.8e %.8e\n", real(S(i,j)), imag(S(i,j)));
       end
       fprintf(fileID, "\n");
   end
   fclose(fileID);
end
