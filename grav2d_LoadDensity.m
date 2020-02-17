%*******************************************************************************
% Function: defRho = grav2d_LoadDensity(pathRho)
% Function: [defRho,defRhoPlane] = grav2d_LoadDensity(pathRho,pathRhoPlane)
%
% Purpose:  Load density variations data in 2D
%
% Inputs:   - pathRho: Path to the depth density variations
%           - pathRhoPlane: Path to the horizontal sectors definition (optional)
%
% Outputs:  - defRho: Structure with the depth-variation density. Two members:
%                     - id: Vector with the density schemes identifiers
%                     - rho: Cell array (same elements as the defRho.id vector)
%                            with the density definitions. Each definition is a
%                            two-columns vector, the first one indicating depth
%                            and the second the density
%           - defRhoPlane: Structure with the horizontal-variation density. Two
%                          members:
%                          - id: Vector with the horizontal schemes identifiers
%                          - zones: Cell array (same elements as the
%                                   defRhoPlane.id vector) with the horizontal
%                                   zones definitions. Each definition is a
%                                   two-column vector, the first one indicating
%                                   the start length and the second the finish
%                                   length
%
% Note: This function does not perform any check about the input arguments
%
% History:  28-06-2018: Function creation
%                       José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************

function [defRho,defRhoPlane] = grav2d_LoadDensity(pathRho,pathRhoPlane)

%Loading density definitions
defRho = aux_LoadData(pathRho);
defRho.rho = defRho.data;
defRho = rmfield(defRho,'data');
%Check if there is horizontal definition
if (nargin==2)&&(nargout>1)
    %Load horizontal definitions
    data_aux = load(pathRhoPlane);
    %Structure creation
    defRhoPlane.id = data_aux(:,1)';
    defRhoPlane.zones = data_aux(:,2:3);
else
    defRhoPlane.id = [];
    defRhoPlane.zones = {};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Auxiliary function for data loading
%Inputs:   - file_path: Working file path
%
%outputs:  - data: Structure with the data. Two members:
%                  - id: Vector with the data identifiers
%                  - data: Cell array (same elements as the data.id vector) with
%                          the correspondent data definitions
function [data] = aux_LoadData(file_path)

%Load data
try
    data_aux = load(file_path);
catch
    error('The file ''%s'' has incorrect format',file_path);
end
%Matrix dimensions
[nr,nc] = size(data_aux);
if (nr<1)||(nc<2)
    error('The file ''%s'' has incorrect format',file_path);
end
%Auxiliary variables
posCont = 1;
posData = 1;
%Infinite loop
while posCont~=0
    %Number of rows and identifier for each definition
    nrd = data_aux(posCont,1);
    id = data_aux(posCont,2);
    %Data extraction
    try
        subdata = data_aux(posCont+1:posCont+nrd,1:2);
    catch
        error('The file ''%s'' has incorrect format',file_path);
    end
    %Assign identifier and data
    data.id(posData) = id;
    data.data{posData} = subdata;
    posData = posData+1;
    %Next header position
    posCont = posCont+nrd+1;
    %Check if the last row was reached
    if posCont>nr
        posCont = 0;
    end
end
%Check if there is any repeated identifier
if length(unique(data.id))<length(data.id)
    error('The file ''%s'' contains repeated identifiers',file_path);
end
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
