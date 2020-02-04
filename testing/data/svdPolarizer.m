function [kl, re, incl, pmod, smod] = svdPolarizer(z, n, e, lambd, SMALL)
    % This is performs the SVD polarization described in Rosenberger's
    % Real-time ground-motion analysis: Distinguishing P and S arrivals
    % in a noisy environment.
    % It takes as input a vertical, north, and east trace (z,n,e) a 
    % decay factor, and a tolerance.  
    % To select a decay factor, lambd, you may
    % try doing the following - Choose a P-wave window, say, 1 s, which at  
    % 100 samples per secondsampling is 100 samples.  Then set lambda 
    % to (Nw-1)/Nw = 99/100.  Likewise, for a 5 second S-wave window
    % you should find have 499/500.  Pretty much - if you make this number
    % greater than or equal to 1 then this algorithm will be metastable
    % or unstable.
    % The tolerance, SMALL, can be made to machine  epsilon if you so 
    % choose but it may be more useful to make close to a fraction of the 
    % standard deviation of the sensor's noise. 
    % The results are the Z, N, E Karhunen-Louve filtered channels,
    % a measure of rectilinearity, a measure of the cosine of the 
    % incidence angle measured from vertical, the P-modulated waveforms
    % and the S-modulated waveforms.  Note, that while the P and S 
    % modulation is applied to the KL filtered data this does not need to
    % be the case and you could instead modulate the original data.
    % Additionally, it may be worthwhile to compute the magnitude of the 
    % modulated signals on top of your original data.
    U = zeros(3,2);
    S = zeros(2,2);
    M = length(z);
    re = zeros(M,1);
    incl = zeros(M,1);
    kl = zeros(M,3);
    % main loop
    for i=1:M
        % Construct data vector
        d = [z(i); n(i); e(i)];
        % Simple initialization
        if (i == 1)
            n2d = norm(d, 2);
            U(:,1) = d/n2d; % Eqn 26
            S(1,1) = n2d;   % Eqn 27
        end
        m = U'*d;  % m || U, Eqn 16
        m_p = d - U*m; % projection _| U, Eq 17
        r = d'*d;
        s = m'*m;
        t = U*m;
        % Innovation p = |(I - U_{n-1} U_{n-1}^T) d_n|
        %              = sqrt( d'd - 2 U'd + d'U U' U U' d )
        %              = sqrt( d'd - 2 U'd + d'U U' d )
        %              = sqrt( r - 2*s + m'm )
        pp = r - 2*s + t'*t;
        abs_pp = abs(pp);
        % Instrument noise level
        if (sqrt(abs_pp) < SMALL)
            abs_pp = 0;
        end
        p = sqrt(abs_pp);
        if (p > 0)
            % Rank increase Q: 3 x 3, Eq. 19
            Q = [lambd*S, m; 0, 0, p];
        else
            % Rank preserving Q: 2 x 3, Eq 22
            Q = [lambd*S, m; 0, 0, 0];
        end
        % Generic SVD for diagonalization Eqn 24
        [A, Ss, ~] = svd(Q);
        % Update S and U
        S = Ss(1:2, 1:2); % Discard smallest singular value
        if (p > 0)
           U = [U, m_p/p]*A(1:3, 1:2); % Rank reduction
        else
           U = U*A(1:2,1:2); % Rank preserving
        end

        % Tabulate rectilinearity, inclination, and KL transform
        re(i) = 0;
        if (S(1,1) > 0)
            re(i) = 1 - S(2,2)/S(1,1); % rectilinearity - modified to match paper
        end
        incl(i) = abs(U(1,1)); % cos(angle of incidence) vs. Z
        kl(i,1:3) = U*U'*d; % kl transform
    end % End main loop
    % Modulate the P and S waves.  There should be a switch in here
    % that requests one or the other.
    pmod = zeros(M, 3);
    smod = zeros(M, 3);
    for j=1:3
        pmod(:,j) = incl.*re.*kl(:,j);
        smod(:,j) = (1 - incl).*re.*kl(:,j);
    end
    
end
    
%SMALL = 0.8;  % 10'th of the standard deviation of the first 200 points
%lambd = 0.99; % 1 second window at 100 samples per second (i.e., 99/100)
%data = load('pb.b206.eh.windowed.txt');
%z = data(:,2);
%n = data(:,3);
%e = data(:,4);

%data = load('pb.b206.3d.txt');    
%shape(data)
%re, incl = svdPolarizer(data(:,2), data(:,3), data(:,4), 1., 1.e-14);
