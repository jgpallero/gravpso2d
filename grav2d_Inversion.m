%PSO INVERSION SCRIPT
clear('all');
%Family members to use (text inside cell): CC, CP, PC, PP, PR, PSO, RC, RP, RR
pso_family_members = {'CP'};
%WORKING FILES
%Prefix for the output files (mandatory)
output_prefix     = 'output_example-nf_';
%Absolute or relative path (folder) for the output files (mandatory)
output_folder     = './examples/output/example-nf/';
%Absolute or relative path to the observations file (mandatory)
observations_file = './examples/data/example-nf/example-nf-profile-utm-observations.txt';
%Absolute or relative path to the initial subsoil file (mandatory)
subsoil_file      = './examples/data/example-nf/example-nf-profile-utm-subsoil.txt';
%Absolute or relative path to the density contrast definition file (mandatory)
density_file      = './examples/data/example-nf/densities.txt';
%Absolute or relative path to the horizontal density contrast file (optional)
h_density_file    = './examples/data/example-nf/density-segments.txt';
%Absolute or relative path to the initial trend file (optional)
trend_file        = './examples/data/example-nf/example-nf-profile-utm-trend.txt';
%Absolute or relative path to the filter window definition file (optional)
filter_file       = './examples/data/example-nf/filter-coefficients.txt';
%Absolute or relative path to the information about boreholes file (optional)
% boreholes_file    = './examples/data/example-nf/boreholes.txt';
%-------------------------------------------------------------------------------
%GENERAL PARAMETERS FOR THE INVERSION
%Use weights in the cost function evaluation: 0/1 -> no/yes (mandatory)
use_weights = 1;
%Swarm size (number of models in each iteration) (mandatory)
swarm_size = 200;
%Iterations number (mandatory)
iterations_number = 150;
%Use only the points on sediments in the inversion: 0/1 -> no/yes (mandatory)
only_points_on_sediments = 1;
%Regional trend estimation in the inversion: 0/1 -> no/yes (mandatory)
%If regional_trend=0 and 'trend_file' is defined, the regional trend is
%substracted from the observations, but their parameters are not optimized
regional_trend = 1;
%Number of times the mean filter is applied (mandatory)
%If the value is <=0 or 'filter_file' is not defined, no filtering is applied
use_filter = 2;
%Size of the filter window (mandatory)
filter_size = 5;
%Weights to the filter based on the prisms' width: 0/1 -> no/yes (mandatory)
%If the prisms have the same width this variable is not used, but its definition
%is mandatory
filter_weight_width = 0;
%Factors (first and second elements) to compute the PSO search space limits, and
%minimum width (third element, in meters) of the search space (mandatory)
%The factors and minimum width are applied to the model stored in
%'subsoil_file', except for the prisms affected by absolute constraints via
%'boreholes_file', which have their depths fixed
depth_factors = [0.50 2.0 10.0];
%Subprism thickness for attraction computation, in meters (optional, 0 default)
%If a value <=0 is assigned, no subprisms are used and the density contrasts
%defined at the surface levels in 'density_file' are used
subprism_size = 20.0;
%Norm to use in the cost function computation (mandatory)
norm_cost_function = 2;
%Minumum and maximum time step. For delta values >1 the algorithm becomes more
%explorative, for values <1 the algorithm becomes more exploitative (mandatory)
%Each row is a couple of minumum (column 1) and maximum (column 2) values
deltats = [0.8 1.2];
%COmputation repetitions for each deltat
repeat = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%THE USER SHOULD NOT MODIFY THE CODE FROM THIS POINT ON%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check deltats definition
if ismatrix(deltats)
    %Dimensions of deltats
    [dtr,dtc] = size(deltats);
    %Check for errors
    if (dtr<1)||(dtc~=2)
        error('The variable ''deltats'' must be a two-column matrix');
    end
else
    error('The variable ''deltats'' must be a matrix');
end
%Check repeat definition
if isscalar(repeat)
    repeat = round(repeat);
    if repeat<1
        error('The variable ''repeat'' must greater or equal than 1');
    end
else
    error('The variable ''repeat'' must be a scalar');
end
%The pso_family_members variable must be a cell array
if iscell(pso_family_members)
    %Loop over all PSO family members
    for i=1:length(pso_family_members)
        %Extract the family member name
        pso_family_member = pso_family_members{i};
        %Prefix for the output file
        prefix_output = sprintf('%s%s_',output_prefix,pso_family_member);
        %Loop over detats values
        for j=1:dtr
            %Deltat values
            deltat = deltats(j,:);
            %Loop over repetitions
            for k=1:repeat
                %Output folder
                folder_output = sprintf('%s/%s/%2.1f-%2.1f-%d/',...
                                        output_folder,pso_family_member,...
                                        deltat(1),deltat(2),k);
                %Print information on screen
                fprintf(1,['\n\nPSO family member: %s, ',...
                           'deltat=[%2.1f %2.1f], repetition=%d\n\n'],...
                        pso_family_member,deltat(1),deltat(2),k);
                %Generate basic structures
                grav2d_BasicStructures;
                %PSO computation
                grav2d_PSO;
            end
        end
    end
else
    error('The variable ''pso_family_members'' must be a cell array');
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
