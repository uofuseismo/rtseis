function envelopeTest()
  t = 0:1/2000:2-1/2000;
  q = chirp(t-2,4,1/2,6,'quadratic',100,'convex').*exp(-4*(t-1).^2);
  q = q + 5; % Add in a bias
  [up,lo] = envelope(q); % Compute envelope with Hilbert transform
  % Write the data
  fileID = fopen('envelopeChirpReference.txt','w');
  for i=1:length(q)
     fprintf(fileID, '%.15e %.13e %.13e\n', q(i), up(i), lo(i));
  end
  end
