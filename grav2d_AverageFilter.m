%*******************************************************************************
% Function: dataf = grav2d_AverageFilter(n,data,cont,fcoef,fixd,fixw,weight)
%           dataf = grav2d_AverageFilter(n,data,cont,fcoef,fixd,fixw,weight,fwc)
%
% Purpose:  Applies a weighted average filter to a data vector
%
% Inputs:   - n: Number of times the filter is applied
%           - data: Vector containing the data to filter
%           - cont: Two-column matrix with the definition of the contiguous data
%                   in the vector 'data':
%                   - Col. 1: Position (referred to 'data') of the first point
%                             for the contiguous set
%                   - Col. 2: Position (referred to 'data') of the last point
%                             for the contiguous set
%           - fcoef: Vector of odd length with the filter window coefficients
%           - fixd: Vector of the same length as 'data' containing identifiers
%                   of fixed data. Two possible values:
%                   - 0: The data is not fixed
%                   - 1: The data is fixed, so the filter coefficient
%                        correspondent to the data will be multiplied by 'fc'
%           - fixw: Weight factor to apply to the filter coefficient
%                   correspondent to the fixed data. Three possible values:
%                   - 0: No special filter is applied to the fixed data
%                   - 1: Previous filtering is applied in the surroundings
%                        (filter window size) of the fixed point prior to the
%                        general filter is applied 'n' times
%                   - Greater than 1:  Previous filtering is applied in the
%                        surroundings (filter window size) of the fixed point
%                        prior to the general filter is applied 'n' times with
%                        the value of 'fixw'
%           - weight: Identifier for weight application:
%                     - 0: No weight application
%                     - Other than 0: Weight application
%           - fwc: Vector of the same length as 'data' containing the weight of
%                  each one. This argument is not mandatory if weight=0
%
% Outputs:  - dataf: Vector containing the filtered data
%
% Note: This function does not perform any check about the input arguments
%
% History:  18-07-2018: Function creation
%                       José Luis García Pallero, jgpallero@gmail.com
%           22-02-2019: Add arguments 'fd' and 'fc', and function reorganization
%                       José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************

function [dataf] = grav2d_AverageFilter(n,data,cont,fcoef,fixd,fixw,weight,fwc)

%Convert data vector to row
data = data(:)';
nd = length(data);
%Convert all vectors to row
fcoef = fcoef(:)';
fixd = fixd(:)';
if weight~=0
    fwc = fwc(:)';
else
    fwc = ones(1,nd);
end
%Number of coefficientes
nc = length(fcoef);
%Output vector
dataf = data;
%Check if the window is 1 length or there is only one data
if (n==0)||(nd<=1)||(nc<=1)
    return;
end
%Half window
hw = floor(nc/2);
%Possible fixed data
pos_fixd = fixd~=0;
%Assign to fixd the fixw value to fixed positions and value 1 to the others
fixd(:) = 1.0;
if fixw~=0.0
    fixd(pos_fixd) = fixw;
end
%Number of subsegments
ns = size(cont,1);
%Loop over subsegments
for i=1:ns
    %Subsegment positions
    pos_s = cont(i,1):cont(i,2);
    nps = length(pos_s);
    %Data correspondent to the subsegment
    data_s = data(pos_s);
    fixd_s = fixd(pos_s);
    pos_fixd_s = pos_fixd(pos_s);
    fwc_s = fwc(pos_s);
    %Check if fixed data exist and prior filtering is needed
    if (fixw>0)&&(sum(pos_fixd_s)>=1)
        %Apply filtering
        ffix = grav2d_AverageFilterAux(data_s,fcoef,fwc_s,fixd_s);
        %Loop over the fixed data
        for j=1:nps
            %If the data is not fixed, next point
            if pos_fixd_s(j)==0
                continue;
            end
            %Loop over half window
            for k=(j-hw):(j+hw)
                %Assign de filtered data
                if (k>=1)&&(k<=nps)
                    if pos_fixd_s(k)==0
                        data_s(k) = ffix(k);
                    end
                end
            end
        end
    end
    %Apply the filter over the data
    for j=1:n
        data_s_f = grav2d_AverageFilterAux(data_s,fcoef,fwc_s,fixd_s);
        %Recover the fixed data
        data_s_f(pos_fixd_s) = data_s(pos_fixd_s);
        %Update the data to filter
        data_s = data_s_f;
    end
    %Assign the result to the output vector
    dataf(pos_s) = data_s;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dfilt] = grav2d_AverageFilterAux(data,fcoef,fweight1,fweight2)

%Output vector
dfilt = data;
%Number of data
nd = length(data);
%Half window
hw = floor(length(fcoef)/2);
%If the filter window is of length 0, the filtered data is the original data
if (nd==1)||(hw==0)
    return;
end
%Check if weights are provided
if nargin<3
    fweight1 = ones(1,nd);
end
if nargin<4
    fweight2 = ones(1,nd);
end
%Zero padding
data = [zeros(1,hw) data zeros(1,hw)];
fweight1 = [zeros(1,hw) fweight1 zeros(1,hw)];
fweight2 = [zeros(1,hw) fweight2 zeros(1,hw)];
%Loop over elements
for i=(hw+1):(hw+nd)
    %Positions in the original vector
    pos = (i-hw):(i+hw);
    %Coefficients
    coef = fcoef.*fweight1(pos).*fweight2(pos);
    %Filtering
    dfilt(i-hw) = sum(coef.*data(pos))/sum(coef);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 2018-2019, J.L.G. Pallero, jgpallero@gmail.com,
%                         J.L. Fernández Martínez, jlfm@uniovi.es
%                         Z. Fernández Muñiz, zulima@uniovi.es
%                         Sylvain Bonvalot, sylvain.bonvalot@ird.fr
%
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%- Redistributions of source code must retain the above copyright notice, this
%  list of conditions and the following disclaimer.
%- Redistributions in binary form must reproduce the above copyright notice,
%  this list of conditions and the following disclaimer in the documentation
%  and/or other materials provided with the distribution.
%- Neither the name of the copyright holders nor the names of its contributors
%  may be used to endorse or promote products derived from this software without
%  specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
%INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
%OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
%ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
