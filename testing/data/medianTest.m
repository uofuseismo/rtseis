function medianTest()
    % Load the reference data
    data = load("gse2.txt");
    n = 11;
    yfilt = medfilt1(data, n);
    % Write the data
    fileID = fopen('medianFilterReference.txt','w');
    % pre-pad to make consistent with IPP
    for i=1:5
        fprintf(fileID, '%.3e\n', 0.0);
    end;
    for i=1:length(yfilt-5)
      fprintf(fileID, '%.3e\n', yfilt(i));
    end;
    fclose(fileID);
end
