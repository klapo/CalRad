% NANEOF_BR3.M
% Function to compute EOFs from demeaned data with NaNs:
% function [B,vars,amps,Nopt] = naneof_br3(data_in);
% Computes EOFs of ROW DIMENSION in d matrix
% are in columns (new!)
% Nopt[1x1] is the number of EOFs generating the least validation error
% EOFs are orthonormal
% variances are in % of total variance ( sum(d(:).^2)/N, not /M/N!)

% following method of Becken and Rixen, JTECH 2003 v20 1839-
%version 2 tries to eliminate column thing
%version 3 optimizes the number of EOFs

%Syntax
% [B,vars,amps,Nopt] = naneof_br3(data_in)
% data_in is of dimension MxN; should have overall nanmean value of 0
% B is of dimension MxM and EOFs are in columns
% vars shows variance described by each EOF and is dimension Mx1
% amps shows the spatial weights, dimension NxM
% Nopt is the number of optimal EOFs for data filling

%Data Reconstruction
% reconstructed = B(:,1:Nopt)*amps(:,1:Nopt)'
% this gives the zero-mean filled data matrix

function [B,vars,amps,Nopt] = naneof_br3(X)


% the data should be de-meaned by now
X0 = X;
idok = find(~isnan(X));

id_val = idok(ceil(rand(floor(length(idok)/50),1)*length(idok))); % validation subset
X_val = X(id_val);
mx_val = sum(sum(X_val.^2));
% remove validation subset from the data
X(id_val) = nan;
idok = find(~isnan(X));


% replace NaNs by zeros
X(isnan(X)) = 0;



Nit = 100;
tol = 1e-5;
% find out how many eigenfunctions to retain...
% dxex = zeros(Nit,10);
% X1 = X;
for Ne=1:min(size(X,2),20)
X1 = X;
    for k=2:Nit
        % compute SVD
        [U,D,V] = svd(X1,0);
        % truncate and estimate "interpolated" D
        %SVD    Singular value decomposition.
%   [U,S,V] = SVD(X) produces a diagonal matrix S, of the same 
%   dimension as X and with nonnegative diagonal elements in
%   decreasing order, and unitary matrices U and V so that
%   X = U*S*V'.
%
%   S = SVD(X) returns a vector containing the singular values.
%
%   [U,S,V] = SVD(X,0) produces the "economy size"
%   decomposition. If X is m-by-n with m > n, then only the
%   first n columns of U are computed and S is n-by-n.
%   For m <= n, SVD(X,0) is equivalent to SVD(X).
%
%   [U,S,V] = SVD(X,'econ') also produces the "economy size"
%   decomposition. If X is m-by-n with m >= n, then it is
%   equivalent to SVD(X,0). For m < n, only the first m columns 
%   of V are computed and S is m-by-m.
%

        N = Ne;
        % truncate
        Ut = U(:,1:N);
        Dt = D(1:N,1:N);
        Vt = V(:,1:N);
        Xa = Ut*Dt*Vt';
        Xa(idok) = X(idok); % restore real data
        X2 = Xa;

        % termination criterium?
        dx=sum((X2-X1).^2,1); %dx = sum(column(X2-X1).^2);
        mx=sum(X2.^2,1); %mx = sum(column(X2.^2));
        dxex = dx/mx;
%         fprintf('Size of dxex is %i .',dxex);
        if dxex <tol,
            fprintf('Converged in %d iterations to the tolerance of %.0e\n',k-1,tol);
            break
        end
        
        X1 = X2;
    end
    % error?
    Xa = Ut*Dt*Vt';
    dx_val = sum(sum((Xa(id_val)-X_val).^2));

    err(Ne) = dx_val/mx_val;
end
% plot(err,'g.-');
% ylabel('Error');
% xlabel('EOFs retained');

Nopt = find(err == min(err), 1);

X1 = X0;
idok = find(~isnan(X1));
X1(isnan(X1)) = 0;
for k=2:Nit
    % compute SVD
    [U,D,V] = svd(X1,0);
    N = Nopt;
    % truncate
    Ut = U(:,1:N);
    Dt = D(1:N,1:N);
    Vt = V(:,1:N);
    Xa = Ut*Dt*Vt';
    Xa(idok) = X(idok); % restore real data
    X2 = Xa;
    
    % termination criterium?
    dx=sum((X2-X1).^2,1); %dx = sum(column(X2-X1).^2);
    mx=sum(X2.^2,1); %mx = sum(column(X2.^2));
    dxex = dx/mx;
    %         fprintf('Size of dxex is %i .',dxex);
    if dxex <tol,
        fprintf('Converged in %d iterations to the tolerance of %.0e\n',k-1,tol);
        break
    end
    
    X1 = X2;
end

% units in B
B = U*D;  %the temporal modes
amps = V; % the spatial modes
% % units in amp
% B = U;
% amps = V*D';

%vars = diag(D)/sum(diag(D))*100;
% total_var=nanmean(X(:).^2); %?
total_var=mean(X(:).^2); %the variance of the original temperature series
vars=diag(D).^2/total_var*100/prod(size(X0)); %divide by the product of the original matrix size and by the variance
