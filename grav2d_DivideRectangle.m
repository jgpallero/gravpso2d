%*******************************************************************************
% Function: division = grav2d_DivideRectangle(rect,rho,st)
%
% Purpose:  Divides a rectangle according its density definition
%
% Inputs:   - rect: Four-element vector containing the geometric definition of
%                   the rectangle (in SI units):
%                   - Pos. 1: Minimum X coordinate
%                   - Pos. 2: Maximum X coordinate
%                   - Pos. 3: Minimum Z coordinate, i.e., Z coordinate of the
%                             rectangle's bottom side
%                   - Pos. 4: Maximum Z coordinate, i.e., Z coordinate of the
%                             rectangle's top side
%           - rho: Two-column matrix containing the density definition:
%                  - Col. 1: Depth, as positive value, from the surface
%                  - Col. 2: Corresponding density for the depth
%           - st: Subrectangle thickness for the original rectangle division. If
%                 a value <=0 is provided only the density at the surface will
%                 be used
%
% Outputs:  - division: Five-column matrix containing the original rectangle
%                       division:
%                       - Col. 1: Minimum X coordinate
%                       - Col. 2: Maximum X coordinate
%                       - Col. 3: Minimum Z coordinate, i.e., Z coordinate of
%                                 the subrectangle's bottom side
%                       - Col. 4: Maximum Z coordinate, i.e., Z coordinate of
%                                 the subrectangle's top side
%                       - Col. 5: Density correspondent to the subrectangle
%
% Note: This function does not perform any check about the input arguments
%
% History:  17-07-2018: Function creation
%                       José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************

function [division] = grav2d_DivideRectangle(rect,rho,st)

%Rows of density definition
nd = size(rho,1);
%Check if the subprism division or the density configuration say that only one
%density must be used
if (nd==1)||(st<=0.0)
    %Rectangle division
    division = [rect rho(1,2)];
else
    %Number of points for interpolation
    npi = 5;
    %Output matrix declaration
    division = [];
    %Depths to absolute heights in rho definition
    rho(:,1) = rect(4)-rho(:,1);
    %Check if there is rectangle above the shallowest rectangle part
    if rect(4)>rho(1,1)
        division = [rect(1:2) rho(1,1) rect(4) rho(1,2)];
    end
    %Loop through the density definition
    for i=1:(nd-1)
        %Block thickness
        tk = rho(i,1)-rho(i+1,1);
        %Rectangle thickness correspondent to the block
        rtk = rho(i,1)-rect(3);
        %Check if we must continue
        if rtk<=0.0
            break;
        elseif rtk>tk
            rtk = tk;
        end
        %Check if the box is constant in density
        if rho(i,2)==rho(i+1,2)
            division = [division;[rect(1:2) rho(i,1)-rtk rho(i,1) rho(i,2)]];
        else
            %Check if the box must be subdivided
            if rtk<st
                %Interpolated density
                irho = interp1([rho(i,1) rho(i+1,1)]',[rho(i,2) rho(i+1,2)]',...
                                linspace(rho(i,1)-rtk,rho(i,1),npi)',...
                                'linear','extrap');
                irho = mean(irho);
                %Add the rectangle
                division = [division;[rect(1:2) rho(i,1)-rtk rho(i,1) irho]];
            else
                %Heights
                hsr = [rho(i,1)-rtk:st:rho(i,1)]';
                hsr = flipud(hsr);
                if hsr(1)<rho(i,1)
                    hsr = [rho(i,1);hsr];
                end
                %Loop over the subrectangles
                for j=1:(length(hsr)-1)
                    %Interpolated density
                    irho = interp1([rho(i,1) rho(i+1,1)]',...
                                   [rho(i,2) rho(i+1,2)]',...
                                   linspace(hsr(j+1),hsr(j),npi)',...
                                   'linear','extrap');
                    irho = mean(irho);
                    %Add the subrectangle
                    division = [division;[rect(1:2) hsr(j+1) hsr(j) irho]];
                end
            end
        end
    end
    %Check if there is rectangle below the deepest rectangle part
    if rect(3)<rho(nd,1)
        division = [division;[rect(1:2) rect(3) rho(nd,1) rho(nd,2)]];
    end
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
