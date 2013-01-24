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
end