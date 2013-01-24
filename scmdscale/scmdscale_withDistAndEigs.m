function [Y, totaltime] = scmdscale_withDistAndEigs(D,p,Ng,Ni,r)
%
%Y = scmdscale(D,p,Ng,Ni,r) takes a data matrix D, and returns an n-by-p 
%   configuration matrix Y.  Rows of Y are the coordinates of n points in
%   p-dimensional space for some p << n. 
%   Ng and Ni are parameters of SC-MDS method.
%   p is an estimated dimension of samples.
%   r = 0: using none random grouping mode; SC-MDS default is radom grouping mode.


	[N m] = size(D);
    

	rand('state',sum(100*clock));
	if r==0
		rand_id = 1:N;
	else
		rand_id = randperm(N);
	end
	
	Nout = p;
    
% 	K = floor((N-Ng)/(Ng-Ni))+1;
		
	X = zeros(Nout,N);
		
	n1 = 1;
	n2 = Ng + mod((N-Ng),(Ng-Ni));
	k = 1;
	tic;
	while (n2 <= N)
		D1 = squareform(pdist(D(rand_id(n1:n2),:)));
        
        X1 = D2X(D1, Nout);
	
		if k == 1
	 		X(:,n1:n2) = X1;
		else
	 		[U b] = moving(X1(:,1:Ni),X(:,n1:(n1+Ni-1)));
% 	 		X1 = U*X1 + kron(b,ones(1,Ng));
            X1 = bsxfun(@plus, U*X1, b);
	 		X(:,n1:n2) = X1;
	 	end
		
		n1 = n2-Ni+1;
		n2 = n2 + Ng - Ni;
		k  = k + 1;
	end
	
	X(:,rand_id) = X;
	X = X';
	X = bsxfun(@minus, X, mean(X,1));
	
	M = X'*X;
	
	[basis,L] = eig(M);
	
	
	Y = X*basis;
	Y = Y(:,end:-1:1);
	totaltime = toc;
return;

function X = C2X(D, n_eig, basics)

%D is correlation matrix

	[m, n] = size(D);
	if m == n
		H = eye(n)-ones(n)/n;
		M = H*D*H;
		[V, DM] = symeig(M, n_eig);

		X = real(diag(sqrt(DM))*V');
	else
		disp('D must be a square matrix');
    end
    
return;    
% end of C2X function

function [X, err] = D2X(D, n_eig, basics)

%D is distance matrix
	[m, n] = size(D);
    
    H = eye(n)-ones(n)/n;
    M = -H*(D.^2)*H/2;
    [V, DM] = symeig(M, n_eig);

    X = real(diag(sqrt(DM))*V');

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
% 		MX  = X - kron(OX, ones(1,nx));
% 		MY  = Y - kron(OY, ones(1,nx));
        MX  = bsxfun(@minus, X, OX);
		MY  = bsxfun(@minus, Y, OY);
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
