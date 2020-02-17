%USER CONFIGURATIONS (all variables are mandatory)
clear('all');
%Absolute or relative path to the observations file
observations_file = './examples/data/example-nf/example-nf-profile-utm.txt';
%Absolute or relative path (folder) for the output files
output_folder = './examples/data/example-nf/';
%Prisms' width identifier for the subsoil model. Three possible values: 0/1/2 ->
%Prisms of variable width/width as mean separation between points/fixed width
idpw = 2;
%Prisms' width (only used if idpw=2)
pw = 500.0;
%Sometimes, in the model definition with idpw=1 or 2 more than one point lies in
%one rectangle. As its top height is fixed to the lowest points, it could be
%some points 'in the free air'. This variable has two possible values:
%- 0: The top prism height is set to the lowest point height, but the other
%     points heights and their gravity anomalies remain unchanged
%- 1: The top prism height is set to the lowest point height, and the other
%     points gravity anomalies are reduced to delete the materials within the
%     point and the top prism's side (the points heights remain unchanged)
reduce_g_top_prisms = 1;
%Initial regional trend estimation. This variable represents the degree of a
%polinomial in the form y=a0+a1*(L-L0)+a2*(L-L0)^2..., where L is the length
%along the profile and L0 a reference length. If rt=0 a constant is estimated,
%rt=1 means a straight line, etc.
%If regional trend estimation is not desired, select rt=-1
rt = 1;
%Use points on sediments for the regional trend estimation: 0/1 -> no/yes
psrt = 0;
%Reduction center for the regional trend, in meters from the FIRST observation
%point (all points included)
%If l0_rt = NaN, it will be computed as the mean of the points lenghts
l0_rt = NaN;
%Density contrast (delta_rho=rho_sediments-rho_basement) for the initial model
%computation. In SI units (kg/m^3)!!!
drho = -590.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%THE USER SHOULD NOT MODIFY THE CODE FROM THIS POINT ON%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK FOR THE EXISTENCE AND KIND OF USER CONFIGURATION VARIABLES
%Observations file
if (exist('observations_file','var')~=1)||(exist(observations_file,'file')~=2)
    error(['The variable ''observations_file'' is not defined or the path ',...
           'provided is not a data file']);
end
%Output folder
if exist('output_folder','var')~=1
    error('The variable ''output_folder'' is not defined');
else
    %Check if the folder exists
    if exist(output_folder,'dir')~=7
        fprintf(2,['Warning: The folder %s does not exist. It will be ',...
                   'created\n'],output_folder);
        mkdir(output_folder);
    end
end
%Prisms' width indentifier
if (exist('idpw','var')~=1)||(isscalar(idpw)==0)
    error('The variable ''idpw'' is not defined or is not a scalar');
end
%Prisms' width
if (idpw==2)&&((exist('pw','var')~=1)||(isscalar(pw)==0))
    error('The variable ''pw'' is not defined or is not a scalar');
end
%Point reduction
if (exist('reduce_g_top_prisms','var')~=1)||(isscalar(reduce_g_top_prisms)==0)
    error(['The variable ''reduce_g_top_prisms'' is not defined or is not a',...
           ' scalar']);
end
%Regional trend
if (exist('rt','var')~=1)||(isscalar(rt)==0)
    error('The variable ''rt'' is not defined or is not a scalar');
end
%Points on sediments
if (exist('psrt','var')~=1)||(isscalar(psrt)==0)
    error('The variable ''psrt'' is not defined or is not a scalar');
end
%Regional trend reduction center
if (exist('l0_rt','var')~=1)||(isscalar(l0_rt)==0)
    error('The variable ''l0_rt'' is not defined or is not a scalar');
end
%Density contrast
if (exist('drho','var')~=1)||(isscalar(drho)==0)
    error('The variable ''drho'' is not defined or is not a scalar');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load observations
obs = load(observations_file);
%Matrix size
[nr,nc] = size(obs);
%Check number of points
if (nr<=2)||(nc<4)
    error('Observations file is not correct (rows<=2 or columns<4)');
end
%The first element of the first row contains the factor to SI units
fsi = obs(1,1);
if fsi==0.0
    error('The conversion factor to SI units must be not zero');
end
%The second element of the first row contains the standard deviations identifier
g_standard_deviations = obs(1,2);
%Check number of columns
if ((g_standard_deviations==0)&&((nc<4)||(nc>5)))||...
   ((g_standard_deviations~=0)&&((nc<5)||(nc>6)))
    error(['Observations file is not correct ',...
           '''g_standard_deviations==0'' and columns<4 or columns>5, or ',...
           '''g_standard_deviations==1'' and columns<5 or columns>6)']);
end
%The variable 'obs' will contain only observations
obs = obs(2:end,:);
nr = nr-1;
%Base name
[obsdir,obsname,obsext] = fileparts(observations_file);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if the observations are given by (x,y) coordinates
if ((g_standard_deviations==0)&&(nc==5))||((g_standard_deviations~=0)&&(nc==6))
    %X coordinates reduced to the centroid
    obs(:,1) = obs(:,1)-mean(obs(:,1));
    %Fit a straight line y=a*x+b to the coordinates
    profile = polyfit(obs(:,1),obs(:,2),1);
    %Points projection onto the line
    u = [1.0 profile(1)];
    u = u./norm(u);
    vecp = obs(:,1:2)-repmat([0.0 profile(2)],nr,1);
    dist = vecp*u';
    %The distances are referred to the first point
    dist = abs(dist-dist(1));
    %The new observations coordinates are the distances
    obs = [dist obs(:,3:nc)];
    nc = nc-1;
else
    %The distances are referred to the first point
    obs(:,1) = obs(:,1)-obs(1,1);
end
%Position of the on-sediments identifier
if g_standard_deviations==0
    psed = 4;
    fprintf(2,['You have indicated the observations file does NOT contain ',...
               'standard deviations\n']);
else
    psed = 5;
    fprintf(2,['You have indicated the observations file CONTAINS standard ',...
               'deviations\n']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check the kind of subsoil partition
if idpw==0
    %Matrix for the subsoil partition
    subsoil = [];
    %Check if the first point is on the sediments
    if obs(1,psed)~=0
        %Distance to the next point
        d = obs(2,1)-obs(1,1);
        %Coordinates for the extremes of the first prism
        subsoil = [obs(1,1)-d/2.0 obs(1,1)+d/2.0 0.0 0.0];
    end
    %Loop through all observed points (except the first and the last one)
    for i=2:(nr-1)
        %Check if the point is on sediments
        if obs(i,psed)~=0
            %Distance to the next point
            d = obs(i+1,1)-obs(i,1);
            %Coordinate for the right extreme of the prism
            subs_aux = [0.0 obs(i,1)+d/2 0.0 0.0];
            %Distance to the previous point (this task were accomplished in the
            %previous step, but in order to maintain clear the code we repeat
            %the task here)
            d = obs(i,1)-obs(i-1,1);
            %Coordinate for the left extreme of the prism
            subs_aux(1) = obs(i,1)-d/2.0;
            %Prism's coordinates
            subsoil = [subsoil;subs_aux];
        end
    end
    %Check if the last point is on the sediments
    if obs(nr,psed)~=0
        %Distance to the previous point
        d = obs(nr,1)-obs(nr-1,1);
        %Coordinates for the extremes of the first prism
        subsoil = [subsoil;[obs(nr,1)-d/2.0 obs(nr,1)+d/2.0 0.0 0.0]];
    end
else
    %Check if the width to use is the mean separation between points
    if idpw==1
        %Prisms mean width
        pw = mean(diff(obs(:,1)));
    end
    %Initial position for the partition
    if obs(1,psed)~=0
        posIni = obs(1,1)-pw/2.0;
    else
        posIni = (obs(1,1)+obs(2,1))/2.0;
    end
    %Final position for the partition
    if obs(nr,psed)~=0
        posFin = obs(nr,1)+pw/2.0;
    else
        posFin = (obs(nr-1,1)+obs(nr,1))/2.0;
    end
    %Number of prisms
    dpos = posFin-posIni;
    np = round(dpos/pw);
    %Correct the initial position
    excess = np*pw-dpos;
    posIni = posIni-excess/2.0;
    %Create the first approximation to the partition
    subsoil = zeros(np,4);
    subsoil(1,1:2) = [posIni posIni+pw];
    for i=2:np
        subsoil(i,1) = subsoil(i-1,2);
        subsoil(i,2) = subsoil(i,1)+pw;
    end
    %Loop for partition adjust
    nIter = 0;
    while 1
        %Iteration number
        nIter = nIter+1;
        %Prisms' centers coordinates
        pcc = (subsoil(:,1)+subsoil(:,2))/2.0;
        %Loop through the observation points
        dm = 0.0;
        for i=1:nr
            %Check if the point is on sediments
            if(obs(i,psed))~=0
                %Distance from the prisms centers to the point
                d = pcc-obs(i,1);
                %Minimum distance in absolute value
                [md,pmd] = min(abs(d));
                %Distance accumulation
                dm = dm+d(pmd);
            end
        end
        %Mean distance between points and prisms' center
        dm = dm/sum(obs(:,psed)~=0);
        %Correct the prisms positions in order to make prisms centers nearer
        %onservation points
        subsoil(:,1:2) = subsoil(:,1:2)-dm;
        %Loop end checking
        if (abs(dm)<1.0e-3)||(nIter==10)
            break;
        end
    end
    %Find points outside sediments
    pos = find(obs(:,psed)==0);
    %The first and last points influence were taken into account previously
    if length(pos)~=0
        if pos(1)==1
            pos = pos(2:end);
        end
        if pos(end)==nr
            pos = pos(1:end-1);
        end
    end
    %Loop through the points
    for i=1:length(pos)
        %Area of influence
        da = obs(pos(i),1)-obs(pos(i)-1,1);
        dp = obs(pos(i)+1,1)-obs(pos(i),1);
        ai = [obs(pos(i),1)-da/2.0 obs(pos(i),1)+dp/2.0];
        %Prisms' centers coordinates
        pcc = (subsoil(:,1)+subsoil(:,2))/2.0;
        %Prisms in the influence area of the point
        pia = (pcc>=ai(1))&(pcc<=ai(2));
        %Delete the affected prisms
        subsoil = subsoil(~pia,:);
    end
end
%Number of prisms
np = size(subsoil,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prisms' centers coordinates
pcc = (subsoil(:,1)+subsoil(:,2))/2.0;
%Height of the upper side of each prism
subsoil(:,3) = interp1(obs(:,1),obs(:,2),pcc,'linear','extrap');
%Reduced points
reduced_points = [];
%Loop through the prisms
for i=1:np
    %Find the points on this prism
    pop = (obs(:,1)>=subsoil(i,1))&(obs(:,1)<=subsoil(i,2));
    npop = sum(pop);
    %The prism height will be the height of the lowest point
    if npop>0
        %Prisms positions
        ppos = find(pop);
        %Top prism height
        h = min(obs(ppos,2));
        %Set the prism height
        subsoil(i,3) = h;
        %Check if the points gravities must be modified
        if (npop>1)&&(reduce_g_top_prisms~=0)
            %points height
            hp = obs(ppos,2);
            for j=1:npop
                %Check if the point is not the lowest
                if (hp(j)-h)~=0.0
                    %Gravity anomaly produced by the subprism below the point
                    dg = grav2d_GravityRectangle([subsoil(i,1:3) hp(j)],drho,...
                                                 obs(ppos(j),1:2))/fsi;
                    %Change the gravity anomaly
                    obs(ppos(j),3) = obs(ppos(j),3)-dg;
                    %Warning message
                    fprintf(2,['Gravity anomaly of point %d was readjusted ',...
                               'in order to take into account the subprism ',...
                               'between it and the prism %d top side\n'],...
                               ppos(j),i);
                    reduced_points = [reduced_points ppos(j)];
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if regional trend must be computed
if rt>=0
    %Reference point for trend computations
    if isnan(l0_rt)
        xc = mean(obs(:,1));
    else
        xc = l0_rt;
    end
    %Coordinates reduced to the centroid
    if psrt~=0
        xrc = obs(:,1)-xc;
        dg = obs(:,3);
    else
        xrc = obs(~logical(obs(:,psed)),1)-xc;
        dg = obs(~logical(obs(:,psed)),3);
    end
    %Check if there are enough points to regional trend determination
    if length(xrc)<rt
        error('There is not enough points to compute the regional trend');
    end
    %Polynomial adjust
    [trend,aux] = polyfit(xrc,dg,rt);
    %Trend values in all the points
    val = polyval(trend,obs(:,1)-xc);
    %Residuals
    res = obs(:,3)-val;
    %Move the regional trend in a such way that the max residual anomaly be 0
    trend(end) = trend(end)+max(res);
    %Trend values in all the points with the displaced parameters
    val = polyval(trend,obs(:,1)-xc);
    %Final residuals
    res = obs(:,3)-val;
    %Check if we are not in GNU Octave or Matlab
    if ~exist('OCTAVE_VERSION')
        %Create the cofactor matrix field in the Matlab structure
        Q = inv(aux.R);
        aux.C = Q*Q';
    end
    %Standard deviations of the trend parameters
    if aux.df~=0
        sd_trend = sqrt(diag(aux.C)/aux.df)*aux.normr;
    else
        sd_trend = zeros(length(trend),1);
    end
else
    %Residuals
    res = obs(:,3);
    %The regional trend values are zero
    val = zeros(nr,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gravity constant (in SI units)
%https://physics.nist.gov/cgi-bin/cuu/Value?bg|search_for=universal_in!
G = 6.67408e-11;
%Residual anomaly for each prism center in SI units
resc = interp1(obs(:,1),res,pcc,'linear','extrap')*fsi;
%Prisms' thickness based on Bouguer plate formula
thickness = resc./(2.0*pi*G*drho);
%Prisms' botton height
subsoil(:,4) = subsoil(:,3)-thickness;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name for the observations file
nof = [output_folder,filesep,obsname,'-observations.txt'];
%Write the observations file
idf = fopen(nof,'wb');
fprintf(idf,'%%Original observations and initial trend substraction\n');
if ~isempty(reduced_points)
    if reduce_g_top_prisms~=0
        fprintf(idf,'%%Gravity anomaly(ies) of point(s)\n%%');
    end
    for i=1:length(reduced_points)
        fprintf(idf,'%d, ',reduced_points(i));
    end
    if reduce_g_top_prisms~=0
        fprintf(idf,'\n%%has(have) been reduced\n');
    end
end
fprintf(idf,'\n');
fprintf(idf,['%%Multiplicative factor to convert gravity into SI units, ',...
             'and standard deviations\n%%identifier (0/1->no/yes). Only ',...
             'the first and second values are used, but the\n%%other ',...
             'values are necessary to match the number of columns\n',...
             '%%       FAC-TO-SI         SD-ID\n']);
if g_standard_deviations==0
    fprintf(idf,['%17.10E             0                  0      0',...
                 '                  0                  0\n'],fsi);
    fprintf(idf,['%%        DISTANCE        HEIGHT            ANOMALY ',...
                 'ON-SED      INITIAL-TREND   RESIDUAL-ANOMALY\n']);
    fprintf(idf,'%17.10E %13.6E %18.11E %6d %18.11E %18.11E\n',[obs val res]');
else
    fprintf(idf,['%17.10E             1                  ',...
                 '0           0      0                  ',...
                 '0                  0\n'],fsi);
    fprintf(idf,['%%        DISTANCE        HEIGHT            ',...
                 'ANOMALY  SD-ANOMALY ON-SED      INITIAL-TREND',...
                 '   RESIDUAL-ANOMALY\n']);
    fprintf(idf,'%17.10E %13.6E %18.11E %11.4E %6d %18.11E %18.11E\n',...
            [obs val res]');
end
fclose(idf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name for the regional trend file
nrtf = [output_folder,filesep,obsname,'-trend.txt'];
%Write the regional trend file
idf = fopen(nrtf,'wb');
if rt<0
    fprintf(idf,'%%There is no regional trend estimation\n');
    fprintf(idf,'%%Do not modify this file\n');
    fprintf(idf,'0\n0\n0\n0\n');
else
    fprintf(idf,['%%Initial regional trend estimation ',...
                 '(residuals=observations-trend)\n']);
    if psrt==0
        fprintf(idf,['%%Only points outside the sediments were used for ',...
                     'the trend computation\n']);
    else
        fprintf(idf,['%%Points on and outside the sediments were used for ',...
                     'the trend computation\n']);
    end
    fprintf(idf,'%%The model is: trend=');
    for i=length(trend):-1:2
        if i>2
            fprintf(idf,'a%d*(L-L0)^%d+',i-1,i-1);
        else
            fprintf(idf,'a1*(L-L0)+');
        end
    end
    fprintf(idf,'a0\n\n');
    fprintf(idf,'%%Multiplicative factor to convert gravity to SI units ');
    fprintf(idf,'(only the first value is used)\n');
    fprintf(idf,'%17.10E',fsi);
    for i=1:length(trend)-1
        fprintf(idf,' 0');
    end
    fprintf(idf,'\n');
    fprintf(idf,'%%Reduction center L0 for the lengths along the profile ');
    fprintf(idf,'(only the first value is used)\n');
    fprintf(idf,'%17.10E',xc);
    for i=1:length(trend)-1
        fprintf(idf,' 0');
    end
    fprintf(idf,'\n');
    fprintf(idf,'%%Regional trend parameters\n%%');
    for i=length(trend):-1:1
        fprintf(idf,' %17s',['a',num2str(i-1)]);
    end
    fprintf(idf,'\n');
    for i=1:length(trend)
        fprintf(idf,' %19.12E',trend(i));
    end
    fprintf(idf,'\n');
    fprintf(idf,'%%Standard deviations of the regional trend parameters\n%%');
    for i=length(trend):-1:1
        fprintf(idf,' %17s',['sd_a',num2str(i-1)]);
    end
    fprintf(idf,'\n');
    for i=1:length(trend)
        fprintf(idf,' %19.12E',sd_trend(i));
    end
    fprintf(idf,'\n');
end
fclose(idf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name for the subsoil file
nsf = [output_folder,filesep,obsname,'-subsoil.txt'];
%Write the observations file
idf = fopen(nsf,'wb');
fprintf(idf,'%%Initial subsoil model computed with delta_rho=%.3f kg/m^3\n',...
        drho);
fprintf(idf,['%%INITIAL-DISTANCE    FINAL-DISTANCE    TOP-HEIGHT ',...
             'BOTTOM-HEIGHT         DEPTH\n']);
fprintf(idf,'%17.10E %17.10E %13.6E %13.6E %13.6E\n',...
        [subsoil subsoil(:,3)-subsoil(:,4)]');
fclose(idf);
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
