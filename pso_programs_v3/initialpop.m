%*******************************************************************************
% Copyright (C), 2018, Juan Luis Fernández Martínez et al., jlfm@uniovi.es
%-------------------------------------------------------------------------------
% Function: [parent] = initialpop(talla,lowlimit,upperlimit)
%
% Purpose:  Generation of random initial population
%
% Inputs:   - talla: Number of models to create
%           - lowlimit: Row vector containing the lower limits for each
%                       parameter
%           - upperlimit: Row vector containing the upper limits for each
%                         parameter
%
% Outputs:  - parent: Matrix of 'talla' rows containing in each row a random
%                     generated model with values between the correspondent
%                     lower and upper limits
%
% Note: This function does not check for the correctness of the input arguments
%-------------------------------------------------------------------------------
%*******************************************************************************

function [parent] = initialpop(talla,lowlimit,upperlimit)

%Number of parameters
nparam = length(lowlimit);
%Range between the limits
rango = upperlimit-lowlimit;
%Matrix of 'talla' rows and a number of columns equal to number of parameters
parent = repmat(lowlimit,talla,1)+repmat(rango,talla,1).*rand(talla,nparam);
