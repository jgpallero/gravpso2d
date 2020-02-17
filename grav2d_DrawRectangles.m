%*******************************************************************************
% Function: h = grav2d_DrawRectangles(rect,col);
%           h = grav2d_DrawRectangles(rect,col,lw);
%
% Purpose:  Draw a set of rectangles and fill them with a solid colour
%
% Inputs:   - rect: Four-column matrix containing the geometric definition of
%                   the rectangles (in SI units):
%                   - Col. 1: Minimum X coordinate
%                   - Col. 2: Maximum X coordinate
%                   - Col. 3: Minimum Z coordinate, i.e., Z coordinate of the
%                             rectangle's bottom side
%                   - Col. 4: Maximum Z coordinate, i.e., Z coordinate of the
%                             rectangle's top side
%           - col: Colour idenificator for each rectangle. Two possibilities:
%                  - Text string with the colour identificator for all the
%                    rectangles or cell array of length equal to the number of
%                    rectangles, containing in each element the text string with
%                    the colour identificator for the correspondent rectangle
%                  - Three-column matrix containing in each one the RGB
%                    component normalized to [0 1] for each rectangle
%           - lw: Line width for the rectangles' border. This argument is
%                 optional
%
% Outputs:  - Plot handle returned by the fill() function used internally
%
% Note: This function does not perform any check about the input arguments
%
% History:  18-08-2018: Function creation
%                       José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************

function [h] = grav2d_DrawRectangles(rect,col,lw)

%Number of rectangles
nr = size(rect,1);
%Colour specification
coltext = 1;
if ischar(col)
    %Create a cell array for the colours
    aux = col;
    col = {};
    for i=1:nr
        col{i} = aux;
    end
elseif isvector(col)
    %Create a matrix
    col = repmat(col,nr,1);
    coltext = 0;
end
%Starting plot order
order = 'h = fill(';
%Loop over rectangles
for i=1:nr
    %Rectangle definition
    aux = sprintf(['[rect(%d,1) rect(%d,2) rect(%d,2) rect(%d,1)],',...
                   '[rect(%d,4) rect(%d,4) rect(%d,3) rect(%d,3)],'],...
                   i,i,i,i,i,i,i,i);
    %Colour definition
    if coltext~=0
        %Colour as text identifier
        aux = sprintf('%s''%s'',',aux,col{i});
    else
        aux = sprintf('%s[%.3f %.3f %.3f],',aux,col(i,1),col(i,2),col(i,3));
    end
    %Add to the global order
    order = [order,aux];
end
%Close the plot order
if nargin>2
    %Add LineWidth
    order = [order,'''LineWidth'',',num2str(lw),');'];
else
    order = order(1:end-1);
    order = [order,');'];
end
%Plot
eval(order);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 2018, J.L.G. Pallero, jgpallero@gmail.com,
%                    J.L. Fernández Martínez, jlfm@uniovi.es
%                    Z. Fernández Muñiz, zulima@uniovi.es
%                    Sylvain Bonvalot, sylvain.bonvalot@ird.fr
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
