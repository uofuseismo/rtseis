function svdPolarizerTest()
   SMALL = 0.8;  % 10'th of the standard deviation of the first 200 points
   lambd = 0.99; % 1 second window at 100 samples per second (i.e., 99/100)
   data = load('pb.b206.eh.windowed.txt');
   z = data(:,2);
   n = data(:,3);
   e = data(:,4);

   [kl, re, incl] = svdPolarizer(z, n, e, lambd, SMALL);
   fileID = fopen('svdReference.txt','w');
   for i=1:length(re)
       fprintf(fileID, '%.8e, %.8e, %.8e, %.13e, %.13e\n', ...
               kl(i,1), kl(i,2), kl(i,3), re(i), incl(i) );
   end
   fclose(fileID);
