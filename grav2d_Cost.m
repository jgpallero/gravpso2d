%*******************************************************************************
% Function: cost = grav2d_Cost(swarm,data,options,opfun)
%
% Purpose:  Computes the cost function in relative percentage
%
% Inputs:   - swarm: Swarm generated by the PSO algorithm
%           - data: 'data' structure
%           - options: 'options' structure
%           - opfun: 'opfun' structure
%
% Outputs:  - cost: Row vector containing the cost function as relative
%                   percentage
%           - grav2d_results: This is a GLOBAL structure, i.e., it is not
%                             assigned explicitly to the function's output
%                             arguments, but must be accessed declaring
%                             'global grav2d_results' in the user script where
%                             this function is called. It contains a field
%                             called 'inv_res' that stores the inversion results
%                             and has this subfields:
%                             - swarm_size: Number of models in each swarm
%                             - iterations: Number of iterations
%                             - obs: Matrix containing the observations
%                                    generated by each model. It has
%                                    swarm_size x iterations rows, and each
%                                    column is a observation point, in the same
%                                    order as it were read in the original data
%                                    file
%                             - f_to_SI_obs: Factor to convert the observations
%                                            into SI units
%                             - model: Matrix of swarm_size x iterations rows
%                                      containing in each one each generated
%                                      model (filtered, if it was indicated)
%                             - trend: Matrix of swarm_size x iterations rows
%                                      containing the regional trend parameters
%                                      If it was no estimated, but used, it
%                                      contains the same parameters repeated for
%                                      each model. If no trend was used, this is
%                                      an empty matrix
%                             - obs_trend: Matrix containing the regional trend
%                                          anomaly in each point. It has
%                                          swarm_size x iterations rows, and
%                                          each column is a observation point,
%                                          in the same order as it were read in
%                                          the original data file
%                             - l0_trend: Base length for the trend computation
%                                         If no trend was used, this is an empty
%                                         matrix
%                             - f_to_SI_trend: Factor to convert the regional
%                                              trend computations into SI units
%                                              If no trend was used, this is an
%                                              empty matrix
%                             - dispersion_f: Vector containing the dispersion
%                                             of the filtered swarm in each
%                                             iteration
%                             - rel_misfit: Vector containing the relative
%                                           misfit, in percentage, for each
%                                           model (containing the regional
%                                           trend, if it was estimated in the
%                                           inversion)
%
% Note: This function does not perform any check about the input arguments
%
% History:  18-07-2018: Function creation
%                       José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************

function [cost] = grav2d_Cost(swarm,data,options,opfun)

%Global variables containing the results
global grav2d_results;
grav2d_results.inv_res.swarm_size = options.pso.size;
grav2d_results.inv_res.iterations = options.pso.maxiter;
%Chek if filtering must be done in the first iteration
if isfield(grav2d_results.extra,'filt_first_it')
    if grav2d_results.extra.filt_first_it~=0
        filter_data = 1;
    else
        filter_data = 0;
    end
else
    filter_data = 1;
end
%Number of models and parameters
[nm,np] = size(swarm);
%Output vector
cost = ones(1,nm);
%Number of prisms
npr = data.subsoil.np;
%Check if regional trend must be taken into account
if (data.trend.use_trend_file~=0)&&(data.trend.regional_trend~=0)
    use_trend = 2;
elseif data.trend.use_trend_file~=0
    use_trend = 1;
    %Trend in SI units
    trend = data.trend.trend0*data.trend.f_to_SI;
else
    use_trend = 0;
    trend = zeros(data.obs.ng,1);
end
%Rectangles matrix
rect = [data.subsoil.l zeros(npr,1) data.subsoil.htop zeros(npr,1)];
rect(:,5) = opfun.subprism_size;
%Filtered models
fmodels = [];
%Initialize observations and observations trend fields in the output structure
if grav2d_results.extra.new_inversion~=0
    grav2d_results.inv_res.obs = [];
    grav2d_results.inv_res.obs_trend = [];
end
%Loop through all models
for i=1:nm
    %Bottom rectangles heights
    brh = swarm(i,1:npr);
    %Chek if filtering must be done in the first iteration
    if filter_data~=0
        %Filtering
        brh = grav2d_AverageFilter(data.filt.use_filter,brh,...
                                   data.subsoil.contiguous,...
                                   data.filt.filter_coef,...
                                   data.subsoil.fixed_depth_f,...
                                   data.subsoil.factor_fixed_depth,...
                                   data.filt.filter_weight_width,...
                                   data.subsoil.filter_weight);
        %Check if any point is below or over the search limits
        over = brh>data.subsoil.hbot_max';
        below = brh<data.subsoil.hbot_min';
        brh(over) = data.subsoil.hbot_max(over)';
        brh(below) = data.subsoil.hbot_min(below)';
    else
        %In the next iteration, models will be filtered (if it was selected)
        grav2d_results.extra.filt_first_it = 1;
    end
    %Add the bottom rectangles heights to the rectangles matrix
    rect(:,3) = brh';
    %Add the heights to the auxiliary filtered models matrix
    fmodels = [fmodels;brh];
    %Gravity generated by the model over all points, in SI units
    grav = grav2d_GravityComputation(rect,data.subsoil.density.rho,data.obs.lh);
    %Regional trend
    if use_trend==2
        %Trend in SI units
        trend = polyval(swarm(i,npr+1:end),data.obs.lh(:,1)-data.trend.l0)*...
                        data.trend.f_to_SI;
    end
    %Substract trend from observations and convert it to working units
    o_t = (data.obs.gSI-trend)*(1.0/data.obs.f_to_SI);
    %Gravity generated by the model in working units
    grav = grav*(1.0/data.obs.f_to_SI);
    %Residuals
    residuals = o_t-grav;
    %Working residuals
    if opfun.only_points_on_sediments~=0
        residuals = residuals(data.obs.ps);
        o_t = o_t(data.obs.ps);
        weights = data.obs.weights(data.obs.ps);
    else
        weights = data.obs.weights;
    end
    %The weights are only applied if the norm is L1 or L2
    if opfun.norm_cost_function==1
        residuals = weights.*residuals;
        o_t = weights.*o_t;
    elseif opfun.norm_cost_function==2
        residuals = sqrt(weights).*residuals;
        o_t = sqrt(weights).*o_t;
    end
    %Cost function
    cost(i) = norm(residuals,opfun.norm_cost_function)/...
              norm(o_t,opfun.norm_cost_function)*100.0;
    %Add the observations to the output structure
    grav2d_results.inv_res.obs = [grav2d_results.inv_res.obs;grav'];
    %Add the observations correspondent to the trend to the output structure
    if use_trend~=0
        grav2d_results.inv_res.obs_trend = [grav2d_results.inv_res.obs_trend
                                            trend'./data.trend.f_to_SI];
    end
end
%Add the factor to convert the observations to SI units
grav2d_results.inv_res.f_to_SI_obs = data.obs.f_to_SI;
%Add the relative misfit to the output structure
if grav2d_results.extra.new_inversion==0
    grav2d_results.inv_res.rel_misfit = [grav2d_results.inv_res.rel_misfit
                                         cost'];
else
    grav2d_results.inv_res.rel_misfit = cost';
end
%Add the models to the output structure
if grav2d_results.extra.new_inversion==0
    grav2d_results.inv_res.model = [grav2d_results.inv_res.model;fmodels];
else
    grav2d_results.inv_res.model = fmodels;
end
%Add the regional trend parameters to the output structure
if use_trend~=0
    if grav2d_results.extra.new_inversion==0
        if use_trend==2
            grav2d_results.inv_res.trend = [grav2d_results.inv_res.trend
                                            swarm(:,npr+1:end)];
            %Add estimated trend parameters to the aux filtered models matrix
            fmodels = [fmodels swarm(:,npr+1:end)];
        else
            grav2d_results.inv_res.trend = [grav2d_results.inv_res.trend
                                            repmat(data.trend.coefficients,...
                                                   nm,1)];
        end
    else
        if use_trend==2
            grav2d_results.inv_res.trend = swarm(:,npr+1:end);
            %Add estimated trend parameters to the aux filtered models matrix
            fmodels = [fmodels swarm(:,npr+1:end)];
        else
            grav2d_results.inv_res.trend = repmat(data.trend.coefficients,nm,1);
        end
    end
    grav2d_results.inv_res.l0_trend = data.trend.l0;
    grav2d_results.inv_res.f_to_SI_trend = data.trend.f_to_SI;
else
    grav2d_results.inv_res.trend = [];
    grav2d_results.inv_res.l0_trend = [];
    grav2d_results.inv_res.f_to_SI_trend = [];
end
%Compute dispersion of filtered models
center_of_gravity = mean(fmodels);
mod_disp2 = (fmodels-repmat(center_of_gravity,nm,1)).^2;
distmed = median(sqrt(sum(mod_disp2')))/norm(center_of_gravity);
%Add dispersion of filtered models
if grav2d_results.extra.new_inversion==0
    grav2d_results.inv_res.dispersion_f = [grav2d_results.inv_res.dispersion_f
                                           distmed];
else
    grav2d_results.inv_res.dispersion_f = distmed;
end
%Check if extra data must be printed
if length(opfun.extra_header)>0
    %Dispersion of filtered models
    fprintf('%8.2f',distmed/grav2d_results.inv_res.dispersion_f(1)*100);
    %Interior models
    modlower = repmat(grav2d_results.model.lowlimit,nm,1);
    modupper = repmat(grav2d_results.model.upperlimit,nm,1);
    minterior = (fmodels>modlower)&(fmodels<modupper);
    minterior = sum(minterior');
    minterior = length(find(minterior>=(np-length(data.filt.filter_coef)+1)));
    minterior = minterior/nm*100;
    fprintf('%11.2f',minterior);
end
%Next time the function is called will not be a new inversion
grav2d_results.extra.new_inversion = 0;
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
