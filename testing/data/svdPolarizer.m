function SVDPolarizer(z, n, e, lambd, SMALL)
    U = zeros(3,2);
    S = zeros(2,2);
    M = length(z);
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
        % Innovation
        pp = r - 2*s + t'*t;
        % Instrument noise level
        if (abs(pp) < SMALL)
            pp = 0;
        end
        p = sqrt(pp);
        if (p > 0)
            % Rank increase Q: 3 x 3, Eq. 19
            Q = [lambd*S, m; 0, 0, p];
        else
            % Rank preserving Q: 2 x 3, Eq 22
            Q = [lambd*S, m;, 0, 0, 0];
        end
        % Generic SVD for diagonalization Eqn 24
        [A, Ss, V] = svd(Q);
        % UPdate S and U
        S = Ss(1:2, 1:2); % Discard smallest singular value
        if (p > 0)
           U = [U, m_p/p]*A(1:3, 1:2); % Rank reduction
        else
           U = U*A(1:2,1:2); % Rank preserving
        end

        re = 1;
        if (S(1,1) > 0)
            re = S(2,2)/S(1,1); % rectilinearit
        end
        incl = acos(abs(U(1,1))); % angle of incidence vs. Z
    end % End main loop
 
