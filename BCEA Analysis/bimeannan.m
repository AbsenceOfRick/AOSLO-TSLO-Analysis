function mu = bimeannan(X,Y,density)
% function mu = bimean(X,Y,density)
% Compute the mean of a bivariate distribution, given the density
% X and Y meshgrid coordinates for the density array
% mu is [E[X], E[Y]]

% 4/2011 bst wrote it

t = nansum(density(:));

mu(1,1) = nansum(nansum(X.*density))/t;
mu(2,1) = nansum(nansum(Y.*density))/t;
