####################################################################################################
# naneof_BR2003.py 
# Karl Lapo July/2015
####################################################################################################
# Function for iteratively calculating eofs on discontinous data 
####################################################################################################
# 	Function to compute EOFs from data with irregular NaNs.
# Computes EOFs of ROW DIMENSION in d matrix are in columns 
# Variances are in # of total variance ( sum(d(:).^2)/N, not /M/N!)
#
# Follows Becken and Rixen, JTECH 2003 http://dx.doi.org/10.1175/1520-0426(2003)020<1839:ECADFF>2.0.CO;2 
#
# SYNTAX:
#	[B,vars,amps,Nopt] = naneof_br3(data_in)
#
# INPUT: 		
#		data_in	= MxN numpy array ----- should have overall nanmean value of 0
#
# OUTPUT:
# 		B 		= MxM numpy array is the PC time series (unitless?), EOF numbers are in columns
# 		vars 	= Mx1 numpy array, shows variance described by each EOF 
# 		amps 	= NXM numpy array, spatial weights (EOFs)
# 		Nopt 	= scalar, the number of optimal EOFs generating the least validation error
#
# Data Reconstruction (this gives the zero-mean filled data matrix):
# 	reconstructed = B(:,1:Nopt)*amps(:,1:Nopt)'

## Import statements
# netcdf/numpy/xray
import numpy as np
from datetime import datetime, timedelta

# OS interaction
import sys, pickle, os

####################################################################################################
# Functions
####################################################################################################

##### Iterative eof 
def naneof_BR2003(X):
	# the data should be de-meaned by now
	X0 = X;
#	idok = find(~isnan(X));
	idok = np.where(not np.isnan(X))

	id_val = idok(ceil(rand(floor(length(idok)/50),1)*length(idok))); # validation subset
	X_val = X(id_val);
	mx_val = sum(sum(X_val.^2));
	# remove validation subset from the data
	X(id_val) = nan;
	idok = find(~isnan(X));


# replace NaNs by zeros
X(isnan(X)) = 0;



Nit = 100;
tol = 1e-5;
# find out how many eigenfunctions to retain...
# dxex = zeros(Nit,10);
# X1 = X;
for Ne=1:min(size(X,2),20)
X1 = X;
    for k=2:Nit
        # compute SVD
        [U,D,V] = svd(X1,0);
        # truncate and estimate "interpolated" D
        #SVD    Singular value decomposition.
#   [U,S,V] = SVD(X) produces a diagonal matrix S, of the same 
#   dimension as X and with nonnegative diagonal elements in
#   decreasing order, and unitary matrices U and V so that
#   X = U*S*V'.
#
#   S = SVD(X) returns a vector containing the singular values.
#
#   [U,S,V] = SVD(X,0) produces the "economy size"
#   decomposition. If X is m-by-n with m > n, then only the
#   first n columns of U are computed and S is n-by-n.
#   For m <= n, SVD(X,0) is equivalent to SVD(X).
#
#   [U,S,V] = SVD(X,'econ') also produces the "economy size"
#   decomposition. If X is m-by-n with m >= n, then it is
#   equivalent to SVD(X,0). For m < n, only the first m columns 
#   of V are computed and S is m-by-m.
#

        N = Ne;
        # truncate
        Ut = U(:,1:N);
        Dt = D(1:N,1:N);
        Vt = V(:,1:N);
        Xa = Ut*Dt*Vt';
        Xa(idok) = X(idok); # restore real data
        X2 = Xa;

        # termination criterium?
        dx=sum((X2-X1).^2,1); #dx = sum(column(X2-X1).^2);
        mx=sum(X2.^2,1); #mx = sum(column(X2.^2));
        dxex = dx/mx;
#         fprintf('Size of dxex is #i .',dxex);
        if dxex <tol,
            fprintf('Converged in #d iterations to the tolerance of #.0e\n',k-1,tol);
            break
        end
        
        X1 = X2;
    end
    # error?
    Xa = Ut*Dt*Vt';
    dx_val = sum(sum((Xa(id_val)-X_val).^2));

    err(Ne) = dx_val/mx_val;
end
# plot(err,'g.-');
# ylabel('Error');
# xlabel('EOFs retained');

Nopt = find(err == min(err), 1);

X1 = X0;
idok = find(~isnan(X1));
X1(isnan(X1)) = 0;
for k=2:Nit
    # compute SVD
    [U,D,V] = svd(X1,0);
    N = Nopt;
    # truncate
    Ut = U(:,1:N);
    Dt = D(1:N,1:N);
    Vt = V(:,1:N);
    Xa = Ut*Dt*Vt';
    Xa(idok) = X(idok); # restore real data
    X2 = Xa;
    
    # termination criterium?
    dx=sum((X2-X1).^2,1); #dx = sum(column(X2-X1).^2);
    mx=sum(X2.^2,1); #mx = sum(column(X2.^2));
    dxex = dx/mx;
    #         fprintf('Size of dxex is #i .',dxex);
    if dxex <tol,
        fprintf('Converged in #d iterations to the tolerance of #.0e\n',k-1,tol);
        break
    end
    
    X1 = X2;
end

# units in B
B = U*D;  #the temporal modes
amps = V; # the spatial modes
# # units in amp
# B = U;
# amps = V*D';

#vars = diag(D)/sum(diag(D))*100;
# total_var=nanmean(X(:).^2); #?
total_var=mean(X(:).^2); #the variance of the original temperature series
vars=diag(D).^2/total_var*100/prod(size(X0)); #divide by the product of the original matrix size and by the variance
