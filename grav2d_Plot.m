%USER CONFIGURATIONS (all variables are mandatory)
clear('all');
%Path (absolute or relative) to the results file (it can be a cell array with
%various files (no more than 20), but in this case they must be referred to the
%same model)
struct_file = {'./examples/output/example-nf/CP/0.8-1.3-1/output_example-nf_CP_grav2d_results.mat',...
               './examples/output/example-nf/CP/0.8-1.3-2/output_example-nf_CP_grav2d_results.mat'};
%Path, absolute or relative, (folder) for the output files
output_folder = './examples/output/example-nf/output-L2/CP/analysis-0.8-1.3';
%Path (absolute or relative) to a file containing a reference model (optional)
ref_model = './examples/data/example-nf/true-model.txt';
%Working equivalent region, in percentage
equivalent_region = 4.5;
%Kind of information to generate: 0/1/2 ->
%only text files/text files and general plots/text files and all plots
kind_of_information = 1;
%Minimum dispersion, in percentage, for swarm collapse
minimum_dispersion = 5.0;
%Plot search bounds: 0/1 -> no/yes
plot_search_bounds = 1;
%Plot format: 'png', 'jpg', 'pdf' or 'svg'
plot_format = 'png';
%-------------------------------------------------------------------------------
%PLOT CONFIGURATIONS
hbin_div_factor = 3;
title_fsize = 14;
xylabel_fsize = 14;
legend_fsize = 11;
ticks_fsize = 11;
line_width = 1.25;
point_size = 6;
obs_size = 20;
obs_symbol = '.';
obs_colour = 'g';
obsn_colour = 'r';
basin_colour = 'b';
equiv_colour = 'm';
bounds_colour = 'r';
ref_colour = 'g';
pos_legend = 'SouthEast'; %best, NorthEast, NorthWest, SouthEast, SouthWest
ltp = {'b^-','rs-','m*-','gd-','co-','bs--','r*--','md--','go--','c^--',...
       'b*-.','rd-.','mo-.','g^-.','cs-.','bd:','ro:','m^:','gs:','c*:'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if the struct_file variable is correct
if (exist('struct_file','var')~=1)||((ischar(struct_file)==0)&&...
   (iscell(struct_file)==0))
    error(['The variable ''struct_file'' is not defined or is not a string ',...
           'nor a cell array']);
end
%Output folders
ofp = [output_folder,filesep,'plots'];
ofpcdf = [ofp,filesep,'cdf'];
oft = [output_folder,filesep,'text'];
ofthist = [oft,filesep,'hist'];
if exist('output_folder','var')~=1
    error('The variable ''output_folder'' is not defined');
else
    %Check if the folder exists
    if exist(output_folder,'dir')~=7
        fprintf(2,['Warning: The folder %s does not exist. It will be ',...
                   'created\n'],output_folder);
        mkdir(output_folder);
        mkdir(oft);
        if kind_of_information~=0
            mkdir(ofthist);
            mkdir(ofp);
            if kind_of_information>1
                mkdir(ofpcdf);
            end
        end
    end
end
%Check if the subfolder 'text' exists
if exist(oft,'dir')~=7
    fprintf(2,['Warning: The folder %s does not exist. It will be ',...
                'created\n'],oft);
    mkdir(oft);
end
if kind_of_information~=0
    %Check if the subfolder 'plots' exists
    if exist(ofp,'dir')~=7
        fprintf(2,['Warning: The folder %s does not exist. It will be ',...
                   'created\n'],ofp);
        mkdir(ofp);
    end
    %Check if the subfolder 'cdf' exists
    if exist(ofthist,'dir')~=7
        fprintf(2,['Warning: The folder %s does not exist. It will be ',...
                   'created\n'],ofthist);
        mkdir(ofthist);
    end
    %Check if the subfolder 'plots/cdf' exists
    if kind_of_information>1
        if exist(ofpcdf,'dir')~=7
            fprintf(2,['Warning: The folder %s does not exist. It will be ',...
                       'created\n'],ofpcdf);
            mkdir(ofpcdf);
        end
    end
end
%Convert the struct_file into a cell array if it is necessary
if ischar(struct_file)
    struct_file = {struct_file};
end
%File loading
nfiles = length(struct_file);
files = {};
for i=1:nfiles
    if exist(struct_file{i},'file')==0
        error('%s is not a file',struct_file{i});
    else
        files{i} = load(struct_file{i});
    end
end
%Reference model
if exist('ref_model','var')==1
    if exist(ref_model,'file')==2
        use_ref_model = 1;
        rmodel = load(ref_model);
        if size(rmodel,2)<2
            error(['The value stored in ''ref_model'' is not a data file\n',...
                   'If you do not want to use a reference model you must ',...
                   'comment the line containing this variable']);
        else
            rmodel = rmodel(:,1:2);
        end
    else
        error(['The value stored in ''ref_model'' is not a data file\n',...
               'If you do not want to use a reference model you must ',...
               'comment the line containing this variable']);
    end
else
    use_ref_model = 0;
    rmodel = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if elements in ltp are more than the number of files to analyze
nltp = length(ltp);
if nfiles>nltp
    %Add elements to ltp repeating the initial series
    for i=nltp+1:nfiles
        ltp{i} = ltp{i-nltp};
    end
end
%Plotting misfit evolution
if kind_of_information~=0
    %This line is here only to force window maximization (works only in Matlab)
    figure('Units','Normalized','Outerposition',[0 0 1 1]);
    clf;
    hold('on');
    legend_text = {};
    max_nit = 0;
    for i=1:nfiles
        nit = files{i}.grav2d_results.inv_res.iterations;
        if nit>max_nit
            max_nit = nit;
        end
        misfit = files{i}.grav2d_results.inv_res.rel_misfit_best';
        [aux,file1,file2] = fileparts(files{i}.grav2d_results.filename);
        if strcmp(upper(files{i}.grav2d_results.options.pso.esquema),'PSO')~=0
            esquema = 'GPSO';
        else
            esquema = [files{i}.grav2d_results.options.pso.esquema,'-PSO'];
        end
        legend_text{i} = [file1,file2,' (',esquema,')'];
        plot_handle = plot(1:nit,misfit,ltp{i});
    end
    hold('off');
    grid('on');
    box('on');
    xlim([0 max_nit+1]);
    set(gca,'FontSize',ticks_fsize); %This before any other FontSize
    title('Convergence curve(s)','FontSize',title_fsize);
    xlabel('Iterations','FontSize',xylabel_fsize);
    ylabel('Best model relative error (%)','FontSize',xylabel_fsize);
    legend(legend_text,'FontSize',legend_fsize,'Interpreter','none');
    set(plot_handle,'LineWidth',line_width);
    set(plot_handle,'MarkerSize',point_size);
    name_file = [ofp,filesep,'convergence.',plot_format];
    if strcmp(plot_format,'pdf')~=0
        orient('Landscape');
        if ~exist('OCTAVE_VERSION')
            print(name_file,['-d',plot_format],'-bestfit');
        else
            print(name_file,['-d',plot_format]);
        end
    else
        print(name_file,['-d',plot_format]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Family identifier
if nfiles>1
    family_id = 'multiple';
else
    if strcmp(upper(files{1}.grav2d_results.options.pso.esquema),'PSO')~=0
        family_id = 'GPSO';
    else
        family_id = [files{1}.grav2d_results.options.pso.esquema,'-PSO'];
    end
end
%Plotting swarm dispersion
if kind_of_information~=0
    clf;
    hold('on');
    for i=1:nfiles
        nit = files{i}.grav2d_results.inv_res.iterations;
        dispersion = files{i}.grav2d_results.inv_res.dispersion_f./...
                     files{i}.grav2d_results.inv_res.dispersion_f(1)*100.0;
        [aux,file1,file2] = fileparts(files{i}.grav2d_results.filename);
        plot_handle = plot([1:nit]',dispersion,ltp{i});
    end
    hold('off');
    grid('on');
    box('on');
    xlim([0 max_nit+1]);
    set(gca,'FontSize',ticks_fsize); %This before any other FontSize
    title('Dispersion curve(s)','FontSize',title_fsize);
    xlabel('Iterations','FontSize',xylabel_fsize);
    ylabel('Filtered swarm relative dispersion (%)','FontSize',xylabel_fsize);
    %legend_text is yet created
    legend(legend_text,'FontSize',legend_fsize,'Interpreter','none');
    set(plot_handle,'LineWidth',line_width);
    set(plot_handle,'MarkerSize',point_size);
    yticks = get(gca,'YTick');
    yticks = [0.0 5.0 10.0 yticks(yticks>10.0)];
    set(gca,'YTick',yticks);
    name_file = [ofp,filesep,'dispersion.',plot_format];
    if strcmp(plot_format,'pdf')~=0
        orient('Landscape');
        if ~exist('OCTAVE_VERSION')
            print(name_file,['-d',plot_format],'-bestfit');
        else
            print(name_file,['-d',plot_format]);
        end
    else
        print(name_file,['-d',plot_format]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search the best model
best_error = Inf;
best_trend = [];
best_model_obst = zeros(1,size(files{1}.grav2d_results.inv_res.obs,2));
for i=1:nfiles
    [mm,pmm] = min(files{i}.grav2d_results.inv_res.rel_misfit);
    if mm<best_error
        best_error = mm;
        best_model = files{i}.grav2d_results.inv_res.model(pmm,:);
        best_model_obs = files{i}.grav2d_results.inv_res.obs(pmm,:);
        best_model_f_SI = files{i}.grav2d_results.inv_res.f_to_SI_obs;
        if ~isempty(files{i}.grav2d_results.inv_res.trend)
            best_trend = files{i}.grav2d_results.inv_res.trend(pmm,:);
            best_trend_l0 = files{i}.grav2d_results.inv_res.l0_trend;
            best_trend_f_SI = files{i}.grav2d_results.inv_res.f_to_SI_trend;
            best_model_obst = files{i}.grav2d_results.inv_res.obs_trend(pmm,:);
        end
    end
end
%Regional trend to file
name_file = [oft,filesep,'best-model-trend.txt'];
idf = fopen(name_file,'wb');
if isempty(best_trend)
    fprintf(idf,'%%There is no regional trend estimation (%s)\n',family_id);
    fprintf(idf,'%%Do not modify this file\n');
    fprintf(idf,'0\n0\n0\n0\n');
else
    nprt = length(best_trend);
    fprintf(idf,'%%Regional trend estimation for the BEST model (%s)\n',...
            family_id);
    fprintf(idf,'%%The model is: trend=');
    for i=nprt:-1:2
        if i>2
            fprintf(idf,'a%d*(L-L0)^%d+',i-1,i-1);
        else
            fprintf(idf,'a1*(L-L0)+');
        end
    end
    fprintf(idf,'a0 (residuals=observations-trend)\n\n');
    fprintf(idf,'%%Multiplicative factor to convert trend to SI units ');
    fprintf(idf,'(only the first value is used)\n');
    fprintf(idf,'%17.10E',best_trend_f_SI);
    for i=1:(nprt-1)
        fprintf(idf,' 0');
    end
    fprintf(idf,'\n');
    fprintf(idf,'%%Reduction center L0 for the lengths along the profile ');
    fprintf(idf,'(only the first value is used)\n');
    fprintf(idf,'%17.10E',best_trend_l0);
    for i=1:(nprt-1)
        fprintf(idf,' 0');
    end
    fprintf(idf,'\n');
    fprintf(idf,'%%Regional trend parameters\n%%');
    for i=nprt:-1:1
        fprintf(idf,' %17s',['a',num2str(i-1)]);
    end
    fprintf(idf,'\n');
    for i=1:nprt
        fprintf(idf,' %19.12E',best_trend(i));
    end
    fprintf(idf,'\n');
    fprintf(idf,'%%No standard deviations are presented in this file\n');
    fprintf(idf,['%%You can determine their uncertainty using the ',...
                 'information contained in the\n',...
                 '%%''*grav2d_results.mat'' file(s)\n%%']);
    for i=nprt:-1:1
        fprintf(idf,' %17s',['sd_a',num2str(i-1)]);
    end
    fprintf(idf,'\n');
    for i=1:nprt
        fprintf(idf,' %19.12E',0.0);
    end
    fprintf(idf,'\n');
end
fclose(idf);
%Observations correspondent to best model to file
name_file = [oft,filesep,'best-model-obs.txt'];
idf = fopen(name_file,'wb');
fprintf(idf,'%%Observations generated by the best model (%s)\n\n',family_id);
fprintf(idf,['%%Multiplicative factor to convert gravity into SI units, ',...
             'and standard deviations\n%%identifier (0/1->no/yes). Only ',...
             'the first and second values are used, but the\n%%other ',...
             'values are necessary to match the number of columns\n',...
             '%%       FAC-TO-SI         SD-ID\n']);
if isempty(files{1}.grav2d_results.data.obs.sd_gSI)
    fprintf(idf,['%17.10E             0                  0      0',...
                 '                  0                  0\n'],best_model_f_SI);
    fprintf(idf,['%%        DISTANCE        HEIGHT            ANOMALY ',...
                 'ON-SED              TREND   RESIDUAL-ANOMALY\n']);
    fprintf(idf,'%17.10E %13.6E %18.11E %6d %18.11E %18.11E\n',...
            [files{1}.grav2d_results.data.obs.lh ...
             best_model_obs'+best_model_obst' ...
             files{1}.grav2d_results.data.obs.ps ...
             best_model_obst' best_model_obs']');
else
    fprintf(idf,['%17.10E             1                  0           0',...
                 '      0                  0                  0\n'],...
            best_model_f_SI);
    fprintf(idf,['%%        DISTANCE        HEIGHT            ANOMALY',...
                 '  SD-ANOMALY ON-SED              TREND',...
                 '   RESIDUAL-ANOMALY\n']);
    fprintf(idf,'%17.10E %13.6E %18.11E %11.4E %6d %18.11E %18.11E\n',...
            [files{1}.grav2d_results.data.obs.lh ...
             best_model_obs'+best_model_obst' ...
             files{1}.grav2d_results.data.obs.sd_gSI(:)./best_model_f_SI ...
             files{1}.grav2d_results.data.obs.ps ...
             best_model_obst' best_model_obs']');
end
fclose(idf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search models inside the working equivalent region
mier = [];
tier = [];
cmier = [];
pos = 0;
for i=1:nfiles
    %Swarm size
    swarm_size = files{i}.grav2d_results.inv_res.swarm_size;
    %Check the existence of valid models
    mft = files{i}.grav2d_results.inv_res.rel_misfit;
    pmft = mft<=equivalent_region;
    mft = mft(pmft);
    if isempty(mft)
        fprintf(2,['Warning: File %s does not contain models inside the ',...
                   'equivalent region\n'],files{i}.grav2d_results.filename);
        return;
    end
    %Loop over iterations
    for j=1:files{i}.grav2d_results.inv_res.iterations
        %Initial position
        pos = (j-1)*swarm_size+1;
        %Misfits
        mft = files{i}.grav2d_results.inv_res.rel_misfit(pos:pos+swarm_size-1);
        %Retain only the misfits inside the working equivalent region
        pmft = mft<=equivalent_region;
        mft = mft(pmft);
        %check if exist models
        if sum(pmft)>0
            %Prisms and trend
            hp = files{i}.grav2d_results.inv_res.model(pos:pos+swarm_size-1,:);
            hp = hp(pmft,:);
            if ~isempty(files{i}.grav2d_results.inv_res.trend)
                tp = ...
                  files{i}.grav2d_results.inv_res.trend(pos:pos+swarm_size-1,:);
                tp = tp(pmft,:);
            else
                tp = [];
            end
            %Check the dispersion
            dispersion_f = files{i}.grav2d_results.inv_res.dispersion_f(j)/...
                           files{i}.grav2d_results.inv_res.dispersion_f(1)*100;
            if dispersion_f<minimum_dispersion
                %Retain the gravity center of the swarm
                if size(hp,1)>1
                    hp = mean(hp);
                    if ~isempty(files{i}.grav2d_results.inv_res.trend)
                        tp = mean(tp);
                    end
                end
            end
            %Add the models to the matrices
            mier = [mier;hp];
            tier = [tier;tp];
            cmier = [cmier;mft];
        end
    end
end
%'Most probable' model and trend
mprob_model = [];
mprob_trend = [];
mprob_limits = [];
%Check if there are working models
[nmier,nprisms] = size(mier);
if nmier==0
    fprintf(2,['Warning: There are no models inside the %.2lf%% equivalent ',...
               'region!'],equivalent_region);
else
    %Prism depths
    pdepth = repmat(files{1}.grav2d_results.data.subsoil.htop',nmier,1)-mier;
    %Check if reference model is used
    if use_ref_model~=0
        %Reference model depths
        hbref = interp1(rmodel(:,1),rmodel(:,2),...
                        files{1}.grav2d_results.data.subsoil.lc,...
                        'linear','extrap');
        depth_ref = files{1}.grav2d_results.data.subsoil.htop-hbref;
    end
    %Numbers of bins for the histograms
    nb = zeros(nprisms,1);
    %Loop over the prisms depths
    for i=1:nprisms
        %Number of histogram bins equals number of meters in depth range divided
        %by hbin_div_factor
        nb(i) = round(abs(max(pdepth(:,i))-min(pdepth(:,i)))/hbin_div_factor);
        if nb(i)<2
            nb(i) = 2;
        end
        %Plot results
        if kind_of_information>1
            clf;
            %Empirical CDF
            subplot(2,1,1);
            [fp,xp] = grav2d_Ecdf(pdepth(:,i));
            hStairs = stairs(xp,fp);
            grid('on');
            set(gca,'FontSize',ticks_fsize); %Before any FontSize configuration
            set(hStairs,'LineWidth',line_width);
            tit = sprintf(['Empirical CDF and histogram, (%d models with ',...
                           'tol<=%.2f%%)\nRectangle %d, L=%.3f m, ',...
                           'top height=%.3f m (%s)'],nmier,equivalent_region,...
                          i,files{1}.grav2d_results.data.subsoil.lc(i),...
                          files{1}.grav2d_results.data.subsoil.htop(i),...
                          family_id);
            title(tit,'FontSize',title_fsize);
            xlabel('Rectangle depth (m)','FontSize',xylabel_fsize);
            ylabel('Probability','FontSize',xylabel_fsize);
            if use_ref_model~=0
                change_lim = 0;
                xl = xlim;
                if (depth_ref(i)<xl(1))||(depth_ref(i)>xl(2))
                    change_lim = 1;
                    xmin = min([xl depth_ref(i)]);
                    xmax = max([xl depth_ref(i)]);
                    inc = (xmax-xmin)/20.0;
                    xl = [min([xl-inc depth_ref(i)-inc]) ...
                          max([xl+inc depth_ref(i)+inc])];
                    xlim(xl);
                end
            end
            %Histogram
            subplot(2,1,2);
            hist(pdepth(:,i),nb(i));
            grid('on');
            set(gca,'FontSize',ticks_fsize); %Before any FontSize configuration
            xlabel('Rectangle depth (m)','FontSize',xylabel_fsize);
            ylabel('Number of occurrences','FontSize',xylabel_fsize);
            %Reference model
            if use_ref_model~=0
                hold('on');
                yl = ylim;
                hRef = plot([depth_ref(i) depth_ref(i)],yl,ref_colour,...
                            'LineWidth',line_width);
                hLeg = legend(hRef,'Reference model','Location',pos_legend);
                if change_lim~=0
                    xlim(xl);
                end
                set(hLeg,'FontSize',legend_fsize);
                hold('off');
            end
            %Print results
            name_file = [ofpcdf,filesep,'cdf-rectangle-',sprintf('%03d',i),...
                         '.',plot_format];
            if strcmp(plot_format,'pdf')~=0
                orient('Portrait');
                if ~exist('OCTAVE_VERSION')
                    print(name_file,['-d',plot_format],'-bestfit');
                else
                    print(name_file,['-d',plot_format]);
                end
            else
                print(name_file,['-d',plot_format]);
            end
        end
        %Height of the most probable bottom part of the prism
        [fh,xh] = hist(pdepth(:,i),nb(i));
        [fhm,posfhm] = max(fh);
        mprob_model = [mprob_model xh(posfhm)];
        %Limits for the 'prob_limits' of models
        prob_limits = 100.0;
        [fp,xp] = grav2d_Ecdf(pdepth(:,i));
        lim1 = interp1(fp,xp,(1.0-prob_limits/100.0)/2.0,'linear','extrap');
        lim2 = interp1(fp,xp,1.0-(1.0-prob_limits/100.0)/2.0,'linear','extrap');
        mprob_limits = [mprob_limits;[lim1 lim2]];
        %Export text file
        name_file = [ofthist,filesep,'hist-rectangle-',sprintf('%03d',i),...
                     '.txt'];
        idf = fopen(name_file,'wb');
        fprintf(idf,'%%Histrogram for rectangle %d (%s)\n',i,family_id);
        fprintf(idf,'%%L=%.3f m, top height=%.3f m\n',...
                files{1}.grav2d_results.data.subsoil.lc(i),...
                files{1}.grav2d_results.data.subsoil.htop(i));
        fprintf(idf,'%%%d models below %.2f%% relative error were used\n\n',...
                length(pdepth(:,i)),equivalent_region);
        fprintf(idf,'%%     DEPTH OCCURRENCES\n');
        fprintf(idf,'%12.6E %10d\n',[xh;fh]);
        fclose(idf);
    end
    mprob_model = files{1}.grav2d_results.data.subsoil.htop'-mprob_model;
    mprob_limits = repmat(files{1}.grav2d_results.data.subsoil.htop,1,2)-...
                   mprob_limits;
    %Check if there are regional trend
    if ~isempty(tier)
        npart = size(tier,2);
        %Number of histogram bins equals the average bins for rectangles
        nb = mean(nb);
        %Loop over the parameters
        for i=1:npart
            %Plot results
            if (kind_of_information>1)&&...
               (files{1}.grav2d_results.data.trend.use_trend_file~=0)&&...
               (files{1}.grav2d_results.data.trend.regional_trend~=0)
                clf;
                %Empirical CDF
                subplot(2,1,1);
                [fh,xh] = grav2d_Ecdf(tier(:,i));
                hStairs = stairs(xh,fh);
                grid('on');
                set(gca,'FontSize',ticks_fsize); %Before any FontSize
                set(hStairs,'LineWidth',line_width);
                tit = sprintf(['Empirical CDF and histogram, (%d models ',...
                               'with tol<=%.2f%%)\nTrend parameter a%d, ',...
                               'L0=%.3f m (%s)'],...
                              size(tier,1),equivalent_region,npart-i,...
                              best_trend_l0,family_id);
                title(tit,'FontSize',title_fsize);
                if (npart-i)==0
                    xltext = sprintf('gravity (m/s^2 * %E)',...
                                  1.0/files{1}.grav2d_results.data.obs.f_to_SI);
                elseif (npart-i)==1
                    xltext = sprintf('gravity/m (gravity=m/s^2 * %E)',...
                                  1.0/files{1}.grav2d_results.data.obs.f_to_SI);
                else
                    xltext = sprintf('gravity/m^%d (gravity=m/s^2 * %E)',...
                                     npart-i,...
                                  1.0/files{1}.grav2d_results.data.obs.f_to_SI);
                end
                xlabel(xltext,'FontSize',xylabel_fsize);
                ylabel('Probability','FontSize',xylabel_fsize);
                %Histogram
                subplot(2,1,2);
                hist(tier(:,i),nb);
                grid('on');
                set(gca,'FontSize',ticks_fsize); %Before any FontSize
                xlabel(xltext,'FontSize',xylabel_fsize);
                ylabel('Number of occurrences','FontSize',xylabel_fsize);
                %Print results
                name_file = [ofpcdf,filesep,'cdf-trend-',...
                             sprintf('%02d',npart-i),'.',plot_format];
                if strcmp(plot_format,'pdf')~=0
                    orient('Portrait');
                    if ~exist('OCTAVE_VERSION')
                        print(name_file,['-d',plot_format],'-bestfit');
                    else
                        print(name_file,['-d',plot_format]);
                    end
                else
                    print(name_file,['-d',plot_format]);
                end
            end
            %Most probable parameters
            [fh,xh] = hist(tier(:,i),nb);
            [fhm,posfhm] = max(fh);
            mprob_trend = [mprob_trend xh(posfhm)];
            %Export text file
            name_file = [ofthist,filesep,'hist-trend-',...
                         sprintf('%02d',npart-i),'.txt'];
            idf = fopen(name_file,'wb');
            fprintf(idf,'%%Histrogram for trend parameter %d (%s)\n',...
                    npart-i,family_id);
            fprintf(idf,'%%%d models below %.2f%% relative error were used\n\n',...
                    length(tier(:,i)),equivalent_region);
            fprintf(idf,'%%     VALUE OCCURRENCES\n');
            fprintf(idf,'%12.6E %10d\n',[xh;fh]);
            fclose(idf);
        end
        %Regional trend to file
        name_file = [oft,filesep,'median-model-trend.txt'];
        idf = fopen(name_file,'wb');
        %Regional trend computation
        trend_mpm = polyval(mprob_trend,...
                            files{1}.grav2d_results.data.obs.lh(:,1)-...
                            best_trend_l0)*...
                            best_trend_f_SI*(1.0/best_model_f_SI);
        nprt = length(mprob_trend);
        fprintf(idf,['%%Regional trend estimation for the MEDIAN model ',...
                     '(%s)\n'],family_id);
        fprintf(idf,'%%The model is: trend=');
        for i=nprt:-1:2
            if i>2
                fprintf(idf,'a%d*(L-L0)^%d+',i-1,i-1);
            else
                fprintf(idf,'a1*(L-L0)+');
            end
        end
        fprintf(idf,'a0 (residuals=observations-trend)\n\n');
        fprintf(idf,'%%Multiplicative factor to convert trend to SI units ');
        fprintf(idf,'(only the first value is used)\n');
        fprintf(idf,'%17.10E',best_trend_f_SI);
        for i=1:(nprt-1)
            fprintf(idf,' 0');
        end
        fprintf(idf,'\n');
        fprintf(idf,'%%Reduction center L0 for the lengths along the profile ');
        fprintf(idf,'(only the first value is used)\n');
        fprintf(idf,'%17.10E',best_trend_l0);
        for i=1:(nprt-1)
            fprintf(idf,' 0');
        end
        fprintf(idf,'\n');
        fprintf(idf,'%%Regional trend parameters\n%%');
        for i=nprt:-1:1
            fprintf(idf,' %17s',['a',num2str(i-1)]);
        end
        fprintf(idf,'\n');
        for i=1:nprt
            fprintf(idf,' %19.12E',mprob_trend(i));
        end
        fprintf(idf,'\n');
        fprintf(idf,'%%No standard deviations are presented in this file\n');
        fprintf(idf,['%%You can determine their uncertainty using the ',...
                    'information contained in the\n',...
                    '%%''*grav2d_results.mat'' file(s)\n%%']);
        for i=nprt:-1:1
            fprintf(idf,' %17s',['sd_a',num2str(i-1)]);
        end
        fprintf(idf,'\n');
        for i=1:nprt
            fprintf(idf,' %19.12E',0.0);
        end
        fprintf(idf,'\n');
    else
        trend_mpm = zeros(size(files{1}.grav2d_results.data.obs.lh,1),1);
        fprintf(idf,'%%There is no regional trend estimation (%s)\n',family_id);
        fprintf(idf,'%%Do not modify this file\n');
        fprintf(idf,'0\n0\n0\n0\n');
    end
    fclose(idf);
    %Gravity generated by the most probable model
    rect = [files{1}.grav2d_results.data.subsoil.l mprob_model' ...
            files{1}.grav2d_results.data.subsoil.htop zeros(nprisms,1)];
    rect(:,5) = files{1}.grav2d_results.opfun.subprism_size;
    grav_mpm = grav2d_GravityComputation(rect,...
                            files{1}.grav2d_results.data.subsoil.density.rho,...
                                         files{1}.grav2d_results.data.obs.lh);
    grav_mpm = grav_mpm*(1.0/best_model_f_SI);
    %Gravity generated by the upper limit
    rect(:,3) = mprob_limits(:,1);
    grav1 = grav2d_GravityComputation(rect,...
                            files{1}.grav2d_results.data.subsoil.density.rho,...
                                      files{1}.grav2d_results.data.obs.lh);
    grav1 = grav1*(1.0/best_model_f_SI);
    %Gravity generated by the lower limit
    rect(:,3) = mprob_limits(:,2);
    grav2 = grav2d_GravityComputation(rect,...
                            files{1}.grav2d_results.data.subsoil.density.rho,...
                                      files{1}.grav2d_results.data.obs.lh);
    grav2 = grav2*(1.0/best_model_f_SI);
    %Gravity generated by the reference model
    if use_ref_model~=0
        rect(:,3) = hbref;
        gravr = grav2d_GravityComputation(rect,...
                            files{1}.grav2d_results.data.subsoil.density.rho,...
                                          files{1}.grav2d_results.data.obs.lh);
        gravr = gravr*(1.0/best_model_f_SI);
    end
    %Substract trend from original observations for best and most probable model
    o_t_best = files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI)-...
               best_model_obst';
    o_t_mpm = files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI)-...
              trend_mpm;
    %Copy of complete vectors
    o_t_best_comp = o_t_best;
    o_t_mpm_comp = o_t_mpm;
    %Residuals
    res_best1 = o_t_best-grav1;
    res_best2 = o_t_best-grav2;
    res_mpm = o_t_mpm-grav_mpm;
    res_mpm1 = o_t_mpm-grav1;
    res_mpm2 = o_t_mpm-grav2;
    res_best_plot = o_t_best-best_model_obs';
    res_mpm_plot = res_mpm;
    if use_ref_model~=0
        res_best_ref = o_t_best-gravr;
        res_mpm_ref = o_t_mpm-gravr;
    end
    %Working residuals
    if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
        pos_ps = files{1}.grav2d_results.data.obs.ps;
        o_t_best = o_t_best(pos_ps);
        o_t_mpm = o_t_mpm(pos_ps);
        res_best1 = res_best1(pos_ps);
        res_best2 = res_best2(pos_ps);
        res_mpm = res_mpm(pos_ps);
        res_mpm1 = res_mpm1(pos_ps);
        res_mpm2 = res_mpm2(pos_ps);
        if use_ref_model~=0
            res_best_ref = res_best_ref(pos_ps);
            res_mpm_ref = res_mpm_ref(pos_ps);
        end
        weights = files{1}.grav2d_results.data.obs.weights(pos_ps);
    else
        weights = files{1}.grav2d_results.data.obs.weights;
    end
    %The weights are only applied if the norm is L1 or L2
    ncf = files{1}.grav2d_results.opfun.norm_cost_function;
    if ncf==1
        res_best1 = weights.*res_best1;
        res_best2 = weights.*res_best2;
        res_mpm = weights.*res_mpm;
        res_mpm1 = weights.*res_mpm1;
        res_mpm2 = weights.*res_mpm2;
        o_t_best = weights.*o_t_best;
        o_t_mpm = weights.*o_t_mpm;
        if use_ref_model~=0
            res_best_ref = weights.*res_best_ref;
            res_mpm_ref = weights.*res_mpm_ref;
        end
    elseif ncf==2
        res_best1 = sqrt(weights).*res_best1;
        res_best2 = sqrt(weights).*res_best2;
        res_mpm = sqrt(weights).*res_mpm;
        res_mpm1 = sqrt(weights).*res_mpm1;
        res_mpm2 = sqrt(weights).*res_mpm2;
        o_t_best = sqrt(weights).*o_t_best;
        o_t_mpm = sqrt(weights).*o_t_mpm;
        if use_ref_model~=0
            res_best_ref = sqrt(weights).*res_best_ref;
            res_mpm_ref = sqrt(weights).*res_mpm_ref;
        end
    end
    %Cost functions
    c_best1 = norm(res_best1,ncf)/norm(o_t_best,ncf)*100.0;
    c_best2 = norm(res_best2,ncf)/norm(o_t_best,ncf)*100.0;
    c_mpm = norm(res_mpm,ncf)/norm(o_t_mpm,ncf)*100.0;
    c_mpm1 = norm(res_mpm1,ncf)/norm(o_t_mpm,ncf)*100.0;
    c_mpm2 = norm(res_mpm2,ncf)/norm(o_t_mpm,ncf)*100.0;
    if use_ref_model~=0
        c_best_ref = norm(res_best_ref,ncf)/norm(o_t_best,ncf)*100.0;
        c_mpm_ref = norm(res_mpm_ref,ncf)/norm(o_t_mpm,ncf)*100.0;
    end
    %Best model to file
    name_file = [oft,filesep,'best-model.txt'];
    idf = fopen(name_file,'wb');
    fprintf(idf,'%%Best model: %.2f%% relative misfit (%s)\n',...
            best_error,family_id);
    fprintf(idf,['%%Statistics generated from all models with relative ',...
                 'misfit below %.2f%%\n'],equivalent_region);
    fprintf(idf,['%%Model composed by all upper limits: %.2f%% relative ',...
                 'misfit\n'],c_best1);
    fprintf(idf,['%%Model composed by all lower limits: %.2f%% relative ',...
                 'misfit\n'],c_best2);
    if files{1}.grav2d_results.opfun.only_points_on_sediments==0
        fprintf(idf,['%%All observation points were involved in the misfit ',...
                     'computation\n']);
    else
        fprintf(idf,['%%Only observation points on the sediments were ',...
                     'involved in the misfit computation\n']);
    end
    fprintf(idf,['%%INITIAL-DISTANCE    FINAL-DISTANCE    TOP-HEIGHT ',...
                 'BOTTOM-HEIGHT   UPPER-LIMIT   LOWER-LIMIT         ',...
                 'DEPTH     MIN-DEPTH     MAX-DEPTH\n']);
    fprintf(idf,['%17.10E %17.10E %13.6E %13.6E %13.6E %13.6E %13.6E ',...
                 '%13.6E %13.6E\n'],...
            [files{1}.grav2d_results.data.subsoil.l ...
             files{1}.grav2d_results.data.subsoil.htop ...
             best_model' mprob_limits ...
             files{1}.grav2d_results.data.subsoil.htop-best_model' ...
             files{1}.grav2d_results.data.subsoil.htop-mprob_limits(:,1) ...
             files{1}.grav2d_results.data.subsoil.htop-mprob_limits(:,2)]');
    fclose(idf);
    %Most probable model to file
    name_file = [oft,filesep,'median-model.txt'];
    idf = fopen(name_file,'wb');
    fprintf(idf,'%%Median model: %.2f%% relative misfit (%s)\n',...
            c_mpm,family_id);
    fprintf(idf,['%%Statistics generated from all models with relative ',...
                 'misfit below %.2f%%\n'],equivalent_region);
    fprintf(idf,['%%Model composed by all upper limits: %.2f%% relative ',...
                 'misfit\n'],c_mpm1);
    fprintf(idf,['%%Model composed by all lower limits: %.2f%% relative ',...
                 'misfit\n'],c_mpm2);
    if files{1}.grav2d_results.opfun.only_points_on_sediments==0
        fprintf(idf,['%%All observation points were involved in the misfit ',...
                     'computation\n']);
    else
        fprintf(idf,['%%Only observation points on the sediments were ',...
                     'involved in the misfit computation\n']);
    end
    fprintf(idf,['%%INITIAL-DISTANCE    FINAL-DISTANCE    TOP-HEIGHT ',...
                 'BOTTOM-HEIGHT   UPPER-LIMIT   LOWER-LIMIT         ',...
                 'DEPTH     MIN-DEPTH     MAX-DEPTH\n']);
    fprintf(idf,['%17.10E %17.10E %13.6E %13.6E %13.6E %13.6E %13.6E ',...
                 '%13.6E %13.6E\n'],...
            [files{1}.grav2d_results.data.subsoil.l ...
             files{1}.grav2d_results.data.subsoil.htop ...
             mprob_model' mprob_limits ...
             files{1}.grav2d_results.data.subsoil.htop-mprob_model' ...
             files{1}.grav2d_results.data.subsoil.htop-mprob_limits(:,1) ...
             files{1}.grav2d_results.data.subsoil.htop-mprob_limits(:,2)]');
    fclose(idf);
    %Observations of most probable model to file
    name_file = [oft,filesep,'median-model-obs.txt'];
    idf = fopen(name_file,'wb');
    fprintf(idf,'%%Observations generated by the median model (%s)\n\n',...
            family_id);
    fprintf(idf,['%%Multiplicative factor to convert gravity into SI ',...
                 'units, and standard deviations\n%%identifier ',...
                 '(0/1->no/yes). Only the first and second values are ',...
                 'used, but the\n%%other values are necessary to match the ',...
                 'number of columns\n',...
                 '%%       FAC-TO-SI         SD-ID\n']);
    if isempty(files{1}.grav2d_results.data.obs.sd_gSI)
        fprintf(idf,['%17.10E             0                  0      0',...
                     '                  0                  0\n'],...
                best_model_f_SI);
        fprintf(idf,['%%        DISTANCE        HEIGHT            ANOMALY ',...
                     'ON-SED              TREND   RESIDUAL-ANOMALY\n']);
        fprintf(idf,'%17.10E %13.6E %18.11E %6d %18.11E %18.11E\n',...
                [files{1}.grav2d_results.data.obs.lh ...
                 grav_mpm+trend_mpm ...
                 files{1}.grav2d_results.data.obs.ps ...
                 trend_mpm grav_mpm]');
    else
        fprintf(idf,['%17.10E             1                  0           0',...
                     '      0                  0                  0\n'],...
                best_model_f_SI);
        fprintf(idf,['%%        DISTANCE        HEIGHT            ANOMALY',...
                     '  SD-ANOMALY ON-SED              TREND',...
                     '   RESIDUAL-ANOMALY\n']);
        fprintf(idf,'%17.10E %13.6E %18.11E %11.4E %6d %18.11E %18.11E\n',...
                [files{1}.grav2d_results.data.obs.lh ...
                 grav_mpm+trend_mpm ...
                 files{1}.grav2d_results.data.obs.sd_gSI(:)./best_model_f_SI ...
                 files{1}.grav2d_results.data.obs.ps ...
                 trend_mpm grav_mpm]');
    end
    fclose(idf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot models and residual anomaly plots
if kind_of_information>0
    %BEST MODEL
    %Minimum and maximum X coordinate, and X limits
    xcoor = [files{1}.grav2d_results.data.obs.lh(:,1)' ...
             files{1}.grav2d_results.data.subsoil.l(1,1) ...
             files{1}.grav2d_results.data.subsoil.l(end,2)];
    xmin = min(xcoor);
    xmax = max(xcoor);
    dx = xmax-xmin;
    dx = dx*0.025;
    %Minimum and maximum Y coordinate, and Y limits
    ycoor = [files{1}.grav2d_results.data.obs.lh(:,2)' mprob_limits(:,2)'];
    if plot_search_bounds~=0
        lowdata = 0;
        upperdata = 0;
        avgt = 0;
        for i=1:nfiles
            avgt1 = mean(files{i}.grav2d_results.model.upperlimit(1:nprisms)-...
                         files{i}.grav2d_results.model.lowlimit(1:nprisms));
            if avgt1>avgt
                avgt = avgt1;
                lowdata = files{i}.grav2d_results.model.lowlimit(1:nprisms);
                upperdata = files{i}.grav2d_results.model.upperlimit(1:nprisms);
            end
        end
        ycoor = [ycoor lowdata];
    end
    ymin = min(ycoor);
    ymax = max(ycoor);
    dy = ymax-ymin;
    dy = dy*0.025;
    %Clear plot window
    clf;
    %Upper window preparation
    p1 = subplot(2,1,1);
    position = get(p1,'Position');
    position(2) = position(2)+0.16;
    position(4) = position(4)-0.15;
    set(p1,'Position',position);
    set(p1,'XTicklabel',[]);
    set(gca,'FontSize',ticks_fsize); %Before any FontSize
    ttitle = sprintf(['Residuals, best model (%.2f%% misfit),\nand %.2f%% ',...
                      'equivalent region of each individual rectangle (%s)'],...
                      best_error,equivalent_region,family_id);
    title(ttitle,'FontSize',title_fsize);
    hold('on');
    box('on');
    grid('on');
    if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
        pos_ps = files{1}.grav2d_results.data.obs.ps;
        hRused = plot(files{1}.grav2d_results.data.obs.lh(pos_ps,1),...
                      res_best_plot(pos_ps),[obs_symbol,obs_colour],...
                      'MarkerSize',obs_size);
        if sum(pos_ps)~=0
            hRnused = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                           res_best_plot(~pos_ps),[obs_symbol,obsn_colour],...
                           'MarkerSize',obs_size);
        end
    else
        hRused = plot(files{1}.grav2d_results.data.obs.lh(:,1),res_best_plot,...
                      [obs_symbol,obs_colour],'MarkerSize',obs_size);
    end
    set(gca,'FontSize',ticks_fsize); %Before any FontSize
    tylabel = sprintf('Residuals\n(m/s^2 * %.1E)',...
                      1.0/files{1}.grav2d_results.data.obs.f_to_SI);
    ylabel(tylabel,'FontSize',xylabel_fsize);
    xlim([xmin-dx xmax+dx]);
    hold('off');
    %Lower window preparation
    p2 = subplot(2,1,2);
    position = get(p2,'Position');
    position(4) = position(4)+0.275;
    set(p2,'Position',position);
    hold('on');
    box('on');
    grid('on');
    %Best model
    rect = [files{1}.grav2d_results.data.subsoil.l ...
            best_model' files{1}.grav2d_results.data.subsoil.htop];
    hBest = grav2d_DrawRectangles(rect,basin_colour,line_width);
    rect = [files{1}.grav2d_results.data.subsoil.l ...
            mprob_limits(:,2) mprob_limits(:,1)];
    hEquiv = grav2d_DrawRectangles(rect,equiv_colour,line_width);
    bmh = [];
    for i=1:nprisms
        aux = [files{1}.grav2d_results.data.subsoil.l(i,1) best_model(i)
               files{1}.grav2d_results.data.subsoil.l(i,2) best_model(i)
               NaN NaN];
        bmh = [bmh;aux];
    end
    hBestL = plot(bmh(:,1),bmh(:,2),basin_colour,'LineWidth',line_width);
    hPoints = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                   files{1}.grav2d_results.data.obs.lh(:,2),...
                   [obs_symbol,obs_colour],'MarkerSize',obs_size);
    if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
        if sum(pos_ps)~=0
            hPointsn = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                            files{1}.grav2d_results.data.obs.lh(~pos_ps,2),...
                            [obs_symbol,obsn_colour],'MarkerSize',obs_size);
        end
    end
    legend_order1 = 'hLeg = legend([hPoints hBest(1) hEquiv(1)';
    legend_order2 = '],''Used points'',''Best model'',''Equivalent region''';
    if plot_search_bounds~=0
        hBounds = plot(files{1}.grav2d_results.data.subsoil.lc,lowdata,...
                       bounds_colour,...
                       files{1}.grav2d_results.data.subsoil.lc,upperdata,...
                       bounds_colour,'LineWidth',line_width);
        hBounds = hBounds(1);
        legend_order1 = [legend_order1,' hBounds'];
        legend_order2 = [legend_order2,',''Search bounds'''];
    else
        hBounds = [];
    end
    if use_ref_model~=0
        hRef = plot(rmodel(:,1),rmodel(:,2),ref_colour,...
                    'LineWidth',line_width+0.5);
        legend_order1 = [legend_order1,' hRef'];
        legend_order2 = [legend_order2,[',''Reference model (',...
                                        sprintf('%.2f%%',c_best_ref),...
                                        ' misfit)''']];
    else
        hRef = [];
    end
    set(gca,'FontSize',ticks_fsize); %Before any FontSize
    ylabel('Height (m)','FontSize',xylabel_fsize);
    xlabel('Length (m)','FontSize',xylabel_fsize);
    xlim([xmin-dx xmax+dx]);
    ylim([ymin-dy ymax+dy]);
    legend_order = [legend_order1,legend_order2,',''Location'',pos_legend);'];
    eval(legend_order);
    set(hLeg,'FontSize',legend_fsize);
    hold('off');
    name_file = [ofp,filesep,'model-best.',plot_format];
    if strcmp(plot_format,'pdf')~=0
        orient('Landscape');
        if ~exist('OCTAVE_VERSION')
            print(name_file,['-d',plot_format],'-bestfit');
        else
            print(name_file,['-d',plot_format]);
        end
    else
        print(name_file,['-d',plot_format]);
    end
    %MOST PROBABLE MODEL
    %Upper window preparation
    p1 = subplot(2,1,1);
    position = get(p1,'Position');
    position(2) = position(2)+0.16;
    position(4) = position(4)-0.15;
    set(p1,'Position',position);
    set(p1,'XTicklabel',[]);
    set(gca,'FontSize',ticks_fsize); %Before any FontSize
    ttitle = sprintf(['Residuals, median model (%.2f%% misfit),\n',...
                      'and %.2f%% equivalent region of each individual ',...
                      'rectangle (%s)'],c_mpm,equivalent_region,family_id);
    title(ttitle,'FontSize',title_fsize);
    hold('on');
    box('on');
    grid('on');
    if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
        pos_ps = files{1}.grav2d_results.data.obs.ps;
        hRused = plot(files{1}.grav2d_results.data.obs.lh(pos_ps,1),...
                      res_mpm_plot(pos_ps),[obs_symbol,obs_colour],...
                      'MarkerSize',obs_size);
        if sum(pos_ps)~=0
            hRnused = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                           res_mpm_plot(~pos_ps),[obs_symbol,obsn_colour],...
                           'MarkerSize',obs_size);
        end
    else
        hRused = plot(files{1}.grav2d_results.data.obs.lh(:,1),res_mpm_plot,...
                      [obs_symbol,obs_colour],'MarkerSize',obs_size);
    end
    set(gca,'FontSize',ticks_fsize); %Before any FontSize
    tylabel = sprintf('Residuals\n(m/s^2 * %.1E)',...
                      1.0/files{1}.grav2d_results.data.obs.f_to_SI);
    ylabel(tylabel,'FontSize',xylabel_fsize);
    xlim([xmin-dx xmax+dx]);
    hold('off');
    %Lower window preparation
    p2 = subplot(2,1,2);
    position = get(p2,'Position');
    position(4) = position(4)+0.275;
    set(p2,'Position',position);
    hold('on');
    box('on');
    grid('on');
    %Best model
    rect = [files{1}.grav2d_results.data.subsoil.l ...
            mprob_model' files{1}.grav2d_results.data.subsoil.htop];
    hBest = grav2d_DrawRectangles(rect,basin_colour,line_width);
    rect = [files{1}.grav2d_results.data.subsoil.l ...
            mprob_limits(:,2) mprob_limits(:,1)];
    hEquiv = grav2d_DrawRectangles(rect,equiv_colour,line_width);
    bmh = [];
    for i=1:nprisms
        aux = [files{1}.grav2d_results.data.subsoil.l(i,1) mprob_model(i)
               files{1}.grav2d_results.data.subsoil.l(i,2) mprob_model(i)
               NaN NaN];
        bmh = [bmh;aux];
    end
    hBestL = plot(bmh(:,1),bmh(:,2),basin_colour,'LineWidth',line_width);
    hPoints = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                   files{1}.grav2d_results.data.obs.lh(:,2),...
                   [obs_symbol,obs_colour],'MarkerSize',obs_size);
    if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
        if sum(pos_ps)~=0
            hPointsn = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                            files{1}.grav2d_results.data.obs.lh(~pos_ps,2),...
                            [obs_symbol,obsn_colour],'MarkerSize',obs_size);
        end
    end
    legend_order1 = 'hLeg = legend([hPoints hBest(1) hEquiv(1)';
    legend_order2 = ['],''Used points'',''Median model'',',...
                     '''Equivalent region'''];
    if plot_search_bounds~=0
        hBounds = plot(files{1}.grav2d_results.data.subsoil.lc,lowdata,...
                       bounds_colour,...
                       files{1}.grav2d_results.data.subsoil.lc,upperdata,...
                       bounds_colour,'LineWidth',line_width);
        hBounds = hBounds(1);
        legend_order1 = [legend_order1,' hBounds'];
        legend_order2 = [legend_order2,',''Search bounds'''];
    else
        hBounds = [];
    end
    if use_ref_model~=0
        hRef = plot(rmodel(:,1),rmodel(:,2),ref_colour,...
                    'LineWidth',line_width+0.5);
        legend_order1 = [legend_order1,' hRef'];
        legend_order2 = [legend_order2,[',''Reference model (',...
                                        sprintf('%.2f%%',c_mpm_ref),...
                                        ' misfit)''']];
    else
        hRef = [];
    end
    set(gca,'FontSize',ticks_fsize); %Before any FontSize
    ylabel('Height (m)','FontSize',xylabel_fsize);
    xlabel('Length (m)','FontSize',xylabel_fsize);
    xlim([xmin-dx xmax+dx]);
    ylim([ymin-dy ymax+dy]);
    legend_order = [legend_order1,legend_order2,',''Location'',pos_legend);'];
    eval(legend_order);
    set(hLeg,'FontSize',legend_fsize);
    hold('off');
    name_file = [ofp,filesep,'model-median.',plot_format];
    if strcmp(plot_format,'pdf')~=0
        orient('Landscape');
        if ~exist('OCTAVE_VERSION')
            print(name_file,['-d',plot_format],'-bestfit');
        else
            print(name_file,['-d',plot_format]);
        end
    else
        print(name_file,['-d',plot_format]);
    end
    %RESIDUAL ANOMALY PLOTS
    %Check if regional trend is used
    if exist('best_model_obst','var')==1
        %Trend values in equal spaced points
        points_trend_p = linspace(files{1}.grav2d_results.data.obs.lh(1,1),...
                                  files{1}.grav2d_results.data.obs.lh(end,1),...
                                  250);
        trend_best_p = polyval(best_trend,points_trend_p-best_trend_l0)*...
                              best_trend_f_SI*(1.0/best_model_f_SI);
        trend_mpm_p = polyval(mprob_trend,points_trend_p-best_trend_l0)*...
                              best_trend_f_SI*(1.0/best_model_f_SI);
        %Margins in Y
        ymax = zeros(1,3);
        ymin = zeros(1,3);
        ymax(1)=max(files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI));
        ymax(2)=max(trend_best_p);
        ymax(3)=max(trend_mpm_p);
        ymin(1)=min(files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI));
        ymin(2)=min(trend_best_p);
        ymin(3)=min(trend_mpm_p);
        ymax = max(ymax);
        ymin = min(ymin);
        dy = (ymax-ymin)*0.025;
        %BEST MODEL
        %Upper window preparation for best model
        p1 = subplot(3,1,1);
        position = get(p1,'Position');
        position(2) = position(2)+0.05;
        position(4) = position(4)-0.05;
        set(p1,'Position',position);
        hold('on');
        box('on');
        grid('on');
        %Plot anomaly
        hAnom = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                  files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI),...
                     basin_colour,'LineWidth',line_width);
        hAnomP = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                  files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI),...
                      [obs_symbol,obs_colour],'MarkerSize',obs_size);
        if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
            if sum(pos_ps)~=0
               hAnomPN = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                              files{1}.grav2d_results.data.obs.gSI(~pos_ps)*...
                              (1.0/best_model_f_SI),...
                              [obs_symbol,obsn_colour],'MarkerSize',obs_size);
            end
        end
        xlim([xmin-dx xmax+dx]);
        ylim([ymin-dy ymax+dy]);
        set(gca,'FontSize',ticks_fsize); %Before any FontSize
        ttitle = sprintf(['Gravity anomalies and trend for the best model ',...
                          '(%s)\nUnits: m/s^2 * %E'],...
                         family_id,1.0/files{1}.grav2d_results.data.obs.f_to_SI);
        title(ttitle,'FontSize',title_fsize);
        ylabel('Original anomaly','FontSize',xylabel_fsize);
        hold('off');
        %Middle window preparation for best model
        p1 = subplot(3,1,2);
        position = get(p1,'Position');
        position(2) = position(2)+0.125;
        position(4) = position(4)-0.05;
        set(p1,'Position',position);
        hold('on');
        box('on');
        grid('on');
        %Plot regional trend
        hAnom = plot(points_trend_p,trend_best_p,basin_colour,...
                     'LineWidth',line_width);
        hAnomP = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                  files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI),...
                      [obs_symbol,obs_colour],...
                      'MarkerSize',obs_size);
        if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
            if sum(pos_ps)~=0
               hAnomPN = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                              best_model_obst(~pos_ps)',...
                              [obs_symbol,obsn_colour],'MarkerSize',obs_size);
            end
        end
        xlim([xmin-dx xmax+dx]);
        ylim([ymin-dy ymax+dy]);
        set(gca,'FontSize',ticks_fsize); %Before any FontSize
        ylabel('Trend','FontSize',xylabel_fsize);
        hold('off');
        %Lower window preparation for best model
        p1 = subplot(3,1,3);
        position = get(p1,'Position');
        position(4) = position(4)+0.15;
        set(p1,'Position',position);
        hold('on');
        box('on');
        grid('on');
        %Plot regional trend
        hAnom = plot(files{1}.grav2d_results.data.obs.lh(:,1),o_t_best_comp,...
                     basin_colour,'LineWidth',line_width);
        hAnomP = plot(files{1}.grav2d_results.data.obs.lh(:,1),o_t_best_comp,...
                      [obs_symbol,obs_colour],'MarkerSize',obs_size);
        if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
            if sum(pos_ps)~=0
               hAnomPN = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                              o_t_best_comp(~pos_ps),...
                              [obs_symbol,obsn_colour],'MarkerSize',obs_size);
            end
        end
        xlim([xmin-dx xmax+dx]);
        set(gca,'FontSize',ticks_fsize); %Before any FontSize
        xlabel('Length (m)','FontSize',xylabel_fsize);
        ylabel('Residual anomaly','FontSize',xylabel_fsize);
        hLeg = legend([hAnomP hAnom],...
                      'Points used in inversion','Gravity anomaly',...
                      'Location',pos_legend);
        set(hLeg,'FontSize',legend_fsize);
        hold('off');
        name_file = [ofp,filesep,'anomaly-best.',plot_format];
        if strcmp(plot_format,'pdf')~=0
            orient('Landscape');
            if ~exist('OCTAVE_VERSION')
                print(name_file,['-d',plot_format],'-bestfit');
            else
                print(name_file,['-d',plot_format]);
            end
        else
            print(name_file,['-d',plot_format]);
        end
        %MOST PROBABLE MODEL
        %Upper window preparation for most probable model
        p1 = subplot(3,1,1);
        position = get(p1,'Position');
        position(2) = position(2)+0.05;
        position(4) = position(4)-0.05;
        set(p1,'Position',position);
        hold('on');
        box('on');
        grid('on');
        %Plot anomaly
        hAnom = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                  files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI),...
                     basin_colour,'LineWidth',line_width);
        hAnomP = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                  files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI),...
                      [obs_symbol,obs_colour],'MarkerSize',obs_size);
        if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
            if sum(pos_ps)~=0
               hAnomPN = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                              files{1}.grav2d_results.data.obs.gSI(~pos_ps)*...
                              (1.0/best_model_f_SI),...
                              [obs_symbol,obsn_colour],'MarkerSize',obs_size);
            end
        end
        xlim([xmin-dx xmax+dx]);
        set(gca,'FontSize',ticks_fsize); %Before any FontSize
        ttitle = sprintf(['Gravity anomalies and trend for the median ',...
                          'model (%s)\nUnits: m/s^2 * %E'],...
                         family_id,1.0/files{1}.grav2d_results.data.obs.f_to_SI);
        title(ttitle,'FontSize',title_fsize);
        ylabel('Original anomaly','FontSize',xylabel_fsize);
        hold('off');
        %Middle window preparation for most probable model
        p1 = subplot(3,1,2);
        position = get(p1,'Position');
        position(2) = position(2)+0.125;
        position(4) = position(4)-0.05;
        set(p1,'Position',position);
        hold('on');
        box('on');
        grid('on');
        %Plot regional trend
        hAnom = plot(points_trend_p,trend_mpm_p,basin_colour,...
                     'LineWidth',line_width);
        hAnomP = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                  files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI),...
                      [obs_symbol,obs_colour],'MarkerSize',obs_size);
        if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
            if sum(pos_ps)~=0
               hAnomPN = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                              trend_mpm(~pos_ps),[obs_symbol,obsn_colour],...
                              'MarkerSize',obs_size);
            end
        end
        xlim([xmin-dx xmax+dx]);
        ylim([ymin-dy ymax+dy]);
        set(gca,'FontSize',ticks_fsize); %Before any FontSize
        ylabel('Trend','FontSize',xylabel_fsize);
        hold('off');
        %Lower window preparation for most probable model
        p1 = subplot(3,1,3);
        position = get(p1,'Position');
        position(4) = position(4)+0.15;
        set(p1,'Position',position);
        hold('on');
        box('on');
        grid('on');
        %Plot regional trend
        hAnom = plot(files{1}.grav2d_results.data.obs.lh(:,1),o_t_mpm_comp,...
                     basin_colour,'LineWidth',line_width);
        hAnomP = plot(files{1}.grav2d_results.data.obs.lh(:,1),o_t_mpm_comp,...
                      [obs_symbol,obs_colour],'MarkerSize',obs_size);
        if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
            if sum(pos_ps)~=0
               hAnomPN = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                              o_t_mpm_comp(~pos_ps),[obs_symbol,obsn_colour],...
                              'MarkerSize',obs_size);
            end
        end
        xlim([xmin-dx xmax+dx]);
        set(gca,'FontSize',ticks_fsize); %Before any FontSize
        xlabel('Length (m)','FontSize',xylabel_fsize);
        ylabel('Residual anomaly','FontSize',xylabel_fsize);
        hLeg = legend([hAnomP hAnom],...
                      'Points used in inversion','Gravity anomaly',...
                      'Location',pos_legend);
        set(hLeg,'FontSize',legend_fsize);
        hold('off');
        name_file = [ofp,filesep,'anomaly-median.',plot_format];
        if strcmp(plot_format,'pdf')~=0
            orient('Landscape');
            if ~exist('OCTAVE_VERSION')
                print(name_file,['-d',plot_format],'-bestfit');
            else
                print(name_file,['-d',plot_format]);
            end
        else
            print(name_file,['-d',plot_format]);
        end
    else
        %If no regional trend is used, the residual anomaly is the used anomaly
        hAnom = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                  files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI),...
                     basin_colour,'LineWidth',line_width);
        hold('on');
        box('on');
        grid('on');
        hAnomP = plot(files{1}.grav2d_results.data.obs.lh(:,1),...
                  files{1}.grav2d_results.data.obs.gSI*(1.0/best_model_f_SI),...
                      [obs_symbol,obs_colour],'MarkerSize',obs_size);
        if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
            if sum(pos_ps)~=0
               hAnomPN = plot(files{1}.grav2d_results.data.obs.lh(~pos_ps,1),...
                              files{1}.grav2d_results.data.obs.gSI(~pos_ps)*...
                              (1.0/best_model_f_SI),...
                              [obs_symbol,obsn_colour],'MarkerSize',obs_size);
            end
        end
        xlim([xmin-dx xmax+dx]);
        set(gca,'FontSize',ticks_fsize); %Before any FontSize
        ttitle = sprintf(['Gravity anomalies for the best and the median ',...
                          'models (%s)\nUnits: m/s^2 * %E'],...
                         family_id,1.0/files{1}.grav2d_results.data.obs.f_to_SI);
        title(ttitle,'FontSize',title_fsize);
        ylabel('Original anomaly','FontSize',xylabel_fsize);
        hold('off');
        set(gca,'FontSize',ticks_fsize); %Before any FontSize
        xlabel('Length (m)','FontSize',xylabel_fsize);
        ylabel('Residual anomaly','FontSize',xylabel_fsize);
        hLeg = legend([hAnomP hAnom],...
                      'Points used in inversion','Gravity anomaly',...
                      'Location',pos_legend);
        set(hLeg,'FontSize',legend_fsize);
        hold('off');
        name_file = [ofp,filesep,'anomaly.',plot_format];
        if strcmp(plot_format,'pdf')~=0
            orient('Landscape');
            if ~exist('OCTAVE_VERSION')
                print(name_file,['-d',plot_format],'-bestfit');
            else
                print(name_file,['-d',plot_format]);
            end
        else
            print(name_file,['-d',plot_format]);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data in the first file
if files{1}.grav2d_results.opfun.only_points_on_sediments~=0
    ps = files{1}.grav2d_results.data.obs.ps;
    working_points1 = files{1}.grav2d_results.data.obs.lh(ps,:);
    working_points_g1 = files{1}.grav2d_results.data.obs.gSI(ps);
else
    working_points1 = files{1}.grav2d_results.data.obs.lh;
    working_points_g1 = files{1}.grav2d_results.data.obs.gSI;
end
subsoil1 = [files{1}.grav2d_results.data.subsoil.l ...
            files{1}.grav2d_results.data.subsoil.htop];
density1 = load(files{1}.grav2d_results.data.density.density_file);
if files{1}.grav2d_results.data.density.use_h_density_file~=0
    hdensity1 = load(files{1}.grav2d_results.data.density.h_density_file);
end
if files{1}.grav2d_results.data.borehole.use_boreholes_file~=0
    borehole1 = load(files{1}.grav2d_results.data.borehole.boreholes_file);
end
%Summary
name_file = [output_folder,filesep,'summary.txt'];
idf = fopen(name_file,'wb');
fprintf(idf,'Files involved in the analysis:\n');
for i=1:nfiles
    fprintf(idf,'%s\n',files{i}.grav2d_results.filename);
end
fprintf(idf,'\n');
for i=1:nfiles
    fprintf(idf,'* %s\n',files{i}.grav2d_results.filename);
    %PSO member, swarm size, and iterations
    fprintf(idf,'  PSO family member:.......................... %s\n',...
            files{i}.grav2d_results.options.pso.esquema);
    fprintf(idf,'  Swarm size:................................. %d\n',...
            files{i}.grav2d_results.options.pso.size);
    fprintf(idf,'  Number of iterations:....................... %d\n',...
            files{i}.grav2d_results.options.pso.maxiter);
    %Norm
    fprintf(idf,'  Norm:....................................... L%d',...
            files{i}.grav2d_results.opfun.norm_cost_function);
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.opfun.norm_cost_function~=...
           files{1}.grav2d_results.opfun.norm_cost_function
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Weights
    fprintf(idf,'  Weights in the cost function computation:... ');
    if files{i}.grav2d_results.data.obs.use_weights~=0
        fprintf(idf,'YES');
    else
        fprintf(idf,'NO');
    end
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.data.obs.use_weights~=...
           files{1}.grav2d_results.data.obs.use_weights
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Used observations
    fprintf(idf,'  Used observations:.......................... ');
    if files{i}.grav2d_results.opfun.only_points_on_sediments==0
        fprintf(idf,'ALL');
    else
        fprintf(idf,'ONLY POINTS ON SEDIMENTS');
    end
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.opfun.only_points_on_sediments~=...
           files{1}.grav2d_results.opfun.only_points_on_sediments
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Regional trend
    fprintf(idf,'  Regional trend:............................. ');
    if files{i}.grav2d_results.data.trend.use_trend_file~=0
        if files{i}.grav2d_results.data.trend.regional_trend~=0
            fprintf(idf,'INVERTED');
        else
            fprintf(idf,'IMPOSED');
        end
        fprintf(idf,', GRADE %d',size(tier,2)-1);
    else
        fprintf(idf,'NO');
    end
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.data.trend.use_trend_file~=...
           files{1}.grav2d_results.data.trend.use_trend_file
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Filtering
    fprintf(idf,'  Filtering:.................................. ');
    if files{i}.grav2d_results.data.filt.use_filter~=0
        fprintf(idf,'YES (%d PASS(ES))',...
                files{i}.grav2d_results.data.filt.use_filter);
    else
        fprintf(idf,'NO');
    end
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.data.filt.use_filter~=...
           files{1}.grav2d_results.data.filt.use_filter
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Filter window width
    fprintf(idf,'  Filter window width:........................ %d',...
            files{i}.grav2d_results.data.filt.filter_size);
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.data.filt.filter_size~=...
           files{1}.grav2d_results.data.filt.filter_size
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Filter weight
    fprintf(idf,'  Filtering weighting using prisms width:..... ');
    if files{i}.grav2d_results.data.filt.filter_weight_width~=0
        fprintf(idf,'YES');
    else
        fprintf(idf,'NO');
    end
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.data.filt.filter_weight_width~=...
           files{1}.grav2d_results.data.filt.filter_weight_width
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Filtering in the first iteration
    fprintf(idf,'  Filtering in the first iteration:........... ');
    if files{i}.grav2d_results.extra.filt_first_it_o~=0
        fprintf(idf,'YES');
    else
        fprintf(idf,'NO');
    end
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.extra.filt_first_it_o~=...
           files{1}.grav2d_results.extra.filt_first_it_o
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Subprisms
    fprintf(idf,'  Subprism size (if applicable):.............. ');
    if files{i}.grav2d_results.opfun.subprism_size==0.0
        fprintf(idf,'NO SUBPRISMS');
    else
        fprintf(idf,'%.2f',files{i}.grav2d_results.opfun.subprism_size);
    end
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.opfun.subprism_size~=...
           files{1}.grav2d_results.opfun.subprism_size
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Horizontal density definition
    fprintf(idf,'  Horizontal density definition:.............. ');
    if files{i}.grav2d_results.data.density.use_h_density_file~=0
        fprintf(idf,'YES');
    else
        fprintf(idf,'NO');
    end
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.data.density.use_h_density_file~=...
           files{1}.grav2d_results.data.density.use_h_density_file
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Boreholes
    fprintf(idf,'  Using boreholes as absolute constraints:.... ');
    if files{i}.grav2d_results.data.borehole.use_boreholes_file~=0
        fprintf(idf,'YES');
    else
        fprintf(idf,'NO');
    end
    if i==1
        fprintf(idf,'\n');
    elseif files{i}.grav2d_results.data.borehole.use_boreholes_file~=0~=...
           files{1}.grav2d_results.data.borehole.use_boreholes_file~=0
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
    else
        fprintf(idf,'\n');
    end
    %Observations file
    fprintf(idf,'  Observations file:.......................... %s',...
            files{i}.grav2d_results.data.obs.observations_file);
    if i==1
        fprintf(idf,'\n');
    else
        if strcmp(files{i}.grav2d_results.data.obs.observations_file,...
                  files{1}.grav2d_results.data.obs.observations_file)==0
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
        end
        if files{i}.grav2d_results.opfun.only_points_on_sediments~=0
            ps = files{i}.grav2d_results.data.obs.ps;
            working_pointsi = files{i}.grav2d_results.data.obs.lh(ps,:);
            working_points_gi = files{i}.grav2d_results.data.obs.gSI(ps);
        else
            working_pointsi = files{i}.grav2d_results.data.obs.lh;
            working_points_gi = files{i}.grav2d_results.data.obs.gSI;
        end
        if (sum(sum(working_pointsi~=working_points1))~=0)||...
           (sum(sum(working_points_gi~=working_points_g1))~=0)
            fprintf(idf,[' *****WARNING: THE POINTS USED ARE NOT THE SAME ',...
                         'AS THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
        else
            fprintf(idf,'\n');
        end
    end
    %Subsoil file
    fprintf(idf,'  Subsoil file:............................... %s',...
            files{i}.grav2d_results.data.subsoil.subsoil_file);
    if i==1
        fprintf(idf,'\n');
    else
        if strcmp(files{i}.grav2d_results.data.subsoil.subsoil_file,...
                  files{1}.grav2d_results.data.subsoil.subsoil_file)==0
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
        end
        subsoili = [files{i}.grav2d_results.data.subsoil.l ...
                    files{i}.grav2d_results.data.subsoil.htop];
        if sum(sum(subsoili~=subsoil1))~=0
            fprintf(idf,[' *****WARNING: THE SUBSOIL DEFINITION USED IS ',...
                         'NOT THE SAME AS THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
        else
            fprintf(idf,'\n');
        end
    end
    %Density file
    fprintf(idf,'  Density definition file:.................... %s',...
            files{i}.grav2d_results.data.density.density_file);
    if i==1
        fprintf(idf,'\n');
    else
        if strcmp(files{i}.grav2d_results.data.density.density_file,...
                  files{1}.grav2d_results.data.density.density_file)==0
            fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE SAME AS ',...
                         'THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
        end
        densityi = load(files{i}.grav2d_results.data.density.density_file);
        if sum(sum(densityi~=density1))~=0
            fprintf(idf,[' *****WARNING: THE DENSITY DEFINITION USED IS ',...
                         'NOT THE SAME AS THE CORRESPONDENT IN %s\n'],...
                    files{1}.grav2d_results.filename);
        else
            fprintf(idf,'\n');
        end
    end
    %Horizontal density file
    if (files{i}.grav2d_results.data.density.use_h_density_file~=0)&&...
       (files{1}.grav2d_results.data.density.use_h_density_file~=0)
        fprintf(idf,'  Horizontal density definition file:......... %s',...
                files{i}.grav2d_results.data.density.h_density_file);
        if i==1
            fprintf(idf,'\n');
        else
            if strcmp(files{i}.grav2d_results.data.density.h_density_file,...
                      files{1}.grav2d_results.data.density.h_density_file)==0
                fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE ',...
                             'SAME AS THE CORRESPONDENT IN %s\n'],...
                        files{1}.grav2d_results.filename);
            end
          hdensityi = load(files{i}.grav2d_results.data.density.h_density_file);
            if sum(sum(hdensityi~=hdensity1))~=0
                fprintf(idf,[' *****WARNING: THE HORIZONTAL DENSITY ',...
                             'DEFINITION USED  IS NOT THE SAME AS THE ',...
                             'CORRESPONDENT IN %s\n'],...
                        files{1}.grav2d_results.filename);
            else
                fprintf(idf,'\n');
            end
        end
    end
    %Filter file
    if (files{i}.grav2d_results.data.filt.use_filter~=0)&&...
       (files{1}.grav2d_results.data.filt.use_filter~=0)
        fprintf(idf,'  Filter file:................................ %s',...
                files{i}.grav2d_results.data.filt.filter_file);
        if i==1
            fprintf(idf,'\n');
        else
            if strcmp(files{i}.grav2d_results.data.filt.filter_file,...
                      files{1}.grav2d_results.data.filt.filter_file)==0
                fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE ',...
                             'SAME AS THE CORRESPONDENT IN %s\n'],...
                        files{1}.grav2d_results.filename);
            end
            if sum(files{i}.grav2d_results.data.filt.filter_coef~=...
                   files{1}.grav2d_results.data.filt.filter_coef)~=0
                fprintf(idf,[' *****WARNING: THE filter DEFINITION USED is ',...
                             'NOT THE SAME AS THE CORRESPONDENT IN %s\n'],...
                        files{1}.grav2d_results.filename);
            else
                fprintf(idf,'\n');
            end
        end
    end
    %Boreholes file
    if (files{i}.grav2d_results.data.borehole.use_boreholes_file~=0)&&...
       (files{1}.grav2d_results.data.borehole.use_boreholes_file~=0)
        fprintf(idf,'  Boreholes file:............................. %s',...
                files{i}.grav2d_results.data.borehole.boreholes_file);
        if i==1
            fprintf(idf,'\n');
        else
            if strcmp(files{i}.grav2d_results.data.borehole.boreholes_file,...
                      files{1}.grav2d_results.data.borehole.boreholes_file)==0
                fprintf(idf,[' *****WARNING: THIS PARAMETER IS NOT THE ',...
                             'SAME AS THE CORRESPONDENT IN %s\n'],...
                        files{1}.grav2d_results.filename);
            end
         boreholei = load(files{i}.grav2d_results.data.borehole.boreholes_file);
            if isequaln(boreholei,borehole1)==0
                fprintf(idf,[' *****WARNING: THE BOREHOLES USED ARE NOT ',...
                             'THE SAME AS THE CORRESPONDENT IN %s\n'],...
                        files{1}.grav2d_results.filename);
            else
                fprintf(idf,'\n');
            end
        end
    end
end
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
