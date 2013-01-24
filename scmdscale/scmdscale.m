function [Y, totaltime] = scmdscale(D,p,Ng,Ni,r)
%
%Y = scmdscale(D,p,Ng,Ni,r) takes an n-by-n distance matrix D, and returns an n-by-p 
%   configuration matrix Y.  Rows of Y are the coordinates of n points in
%   p-dimensional space for some p << n. 
%   Ng and Ni are parameters of SC-MDS method.
%   p is an estimated dimension of samples.
%   r = 0: using none random grouping mode; SC-MDS default is radom grouping mode.


	[n m] = size(D);
	del = 10*eps(class(D));
	
	% lower triangle form for D, make sure it's a valid dissimilarity matrix
	if n == 1
	    n = (1+sqrt(1+8*m))/2;
	    if n == fix(n) && all(D >= 0)
	        D = squareform(D); % assumes zero diagonal, similarity not allowed
	        [n m] = size(D);
	        
	        if all(diag(D)<del)
	        	matrix_mode = 'diss';
	        elseif all(abs(diag(D) - 1) < del) && all(all(D < 1+del))
	        	matrix_mode = 'corr';
            else
                error('scmdscale:BadDistance',...
                  'Not a valid dissimilarity or similarity matrix.')
            end
	
        end
	    
	% full matrix form, make sure it's valid similarity/dissimilarity matrix
	elseif n == m && all(all(D >= 0 & abs(D - D') <= del*max(max(D))))
	    
	    % it's a dissimilarity matrix
	    if all(diag(D) < del)
	        % nothing to do
	        matrix_mode = 'diss';
	    % it's a similarity matrix -- transform to dissimilarity matrix.
	    % the sqrt is not entirely arbitrary, see Seber, eqn. 5.73
	    elseif all(abs(diag(D) - 1) < del) && all(all(D < 1+del))
	        matrix_mode = 'diss';
	        D = 1 - D;
	    else
	        matrix_mode = 'tenser_product';    
	    end
	else
	    error('scmdscale:BadDistance',...
	          'Not a valid dissimilarity or distance matrix.')
	end
	
	N = n;
	if nargin == 4
		r = 1;
	end
	
	if nargin ==2
		r = 1;
		Ni = p+1;
		Ng = 2*Ni;
	end
	
	if nargin < 2
		error('scmdscale: no enough input data.')
	end
	
	rand('state',sum(100*clock));
	if r==0
		rand_id = 1:N;
	else
		rand_id = randperm(N);
	end
	
	Nout = p;
	K = floor((N-Ng)/(Ng-Ni))+1;
		
	X = zeros(Nout,N);
		
	n1 = 1;
	n2 = Ng + mod((N-Ng),(Ng-Ni));
	k = 1;
	tic;
	while (n2 <= N)
		D1 = D(rand_id(n1:n2),rand_id(n1:n2));
		if strcmp(matrix_mode,'diss')
			X1 = D2X(D1, Nout);
		else
			X1 = C2X(D1,Nout);
		end
	
		if k == 1
	 		X(:,n1:n2) = X1;
		else
	 		[U b] = moving(X1(:,1:Ni),X(:,n1:(n1+Ni-1)));
	 		X1 = U*X1 + kron(b,ones(1,Ng));
	 		X(:,n1:n2) = X1;
	 	end
		
		n1 = n2-Ni+1;
		n2 = n2 + Ng - Ni;
		k  = k + 1;
	end
	
	X(:,rand_id) = X;
	X = X';
	X = zero_sum(X);
	
	M = X'*X;
	
	[basis,L] = eig(M);
	
	
	Y = X*basis;
	Y = Y(:,end:-1:1);
	totaltime = toc;
return;

function X = C2X(D, n_eig, basics)

%D is correlation matrix

	[m, n] = size(D);
	TOL = 1e-6;
	if m == n
		H = eye(n)-ones(n)/n;
		M = H*D*H;
		[V, DM] = eig(M);
		DM = real(diag(DM));
		if nargin == 1
			ID = find(DM > TOL);
            sum(DM>TOL)
		elseif nargin == 2
			[Y ID] = sort(DM,1,'descend');
			ID = ID(1:n_eig);
		elseif nargin == 3
			[p q] = size(basics);
			[Y ID] = sort(DM,1,'descend');
			ID = ID(1:p);			
		end
		V = V(:,ID);
		DM = DM(ID);
		X = real(sqrt(diag(DM))*V');
		if nargin == 3
			X = basics'*X;
		end
	else
		disp('D must be a square matrix');
    end
    
return;    
% end of C2X function

function [X, err] = D2X(D, n_eig, basics)

%D is distance matrix

	[m, n] = size(D);
	TOL = 1e-8;
	if m == n
		H = eye(n)-ones(n)/n;
		M = -H*(D.^2)*H/2;
		[V, DM] = eig(M);
		DM = real(diag(DM));
		if nargin == 1
			ID = find(DM > TOL);
            sum(DM>TOL)
		elseif nargin == 2
			[Y ID] = sort(DM,1,'descend');
			ID = ID(1:n_eig);
			if n_eig == length(Y)
				err = 0;
			else
				err = abs(Y(n_eig+1));
            end
		elseif nargin == 3
			[p q] = size(basics);
			[Y ID] = sort(DM,1,'descend');
			ID = ID(1:p);			
		end
		V = V(:,ID);
		DM = DM(ID);
		X = real(sqrt(diag(DM))*V');
		if nargin == 3
			X = basics'*X;
		end
	else
		disp('D must be a square matrix');
    end

return;
% end of D2X function;

function [U b] = moving(X, Y)
% Return a affine map Y = UX + kron(b,ones(1,n)), UU' = I
% X = [x_1, x_2, ... , x_n]; x_j \in R^m
% Y = [y_1, y_2, ... , y_n]; y_j \in R^m

	[mx nx] = size(X);
	[my ny] = size(Y);
	if mx == my && nx == ny
		%OX = X(:,1);
		%OY = Y(:,1);
		OX = sum(X,2)/nx;
		OY = sum(Y,2)/ny;
		MX  = X - kron(OX, ones(1,nx));
		MY  = Y - kron(OY, ones(1,nx));
		[QX RX] = qr(MX,0); % QX' * (X-OX) = R = QY' * (Y-OY)
		[QY RY] = qr(MY,0); % U = QY*QX'        
		[m, n]  = size(QX);
		for p=1:n
            if sign(RX(p,p))~=sign(RY(p,p))
			    QY(:,p) = -QY(:,p);
			    %RY(p,:) = -RY(p,:); 
            end
		end
		U = QY*QX';
		b = OY-U*OX;
		%max(max(abs(Y-U*X-kron(b,ones(1,nx)))))
	else
		disp('Matrix size of X is not equal to Y');
		U = zeros(mx,mx);
		b = zeros(mx,1);
	end
return;
%end of moving function;

function M = zero_sum(A)

	[m, n] = size(A);
	M = A-repmat(mean(A,1),m,1);
return;
%end of zero_sum function