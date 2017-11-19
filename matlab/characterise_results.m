function characterise_results(fname)
%CHARACTERISE_RESULTS Plots characteristics of COMSOL results `fname`
%
%DESCRIPTION
%
%
%INPUTS
%   fname - COMSOL results .txt file to parse
%
%COPYRIGHT (C) Russell Maguire 2017
[paramv, Tv, C, dz, cte, paramname] = parse_fname(fname);
Tname = 'Temperature (\si{\celsius})';
Cname = 'Capacitance (\si{\femto\farad})';
dzname = 'Displacement (\si{\micro\meter})';
ctename = 'Coefficient of thermal expansion (\si{\per\kelvin})';
stressname = '\ce{SiO2} pre-stress (\si\{\giga\pascal})';
if strcmp(fname, 'material_study.txt')
    paramname = ctename;
    paramv = cte(1,:);
elseif strcmp(fname, 'prestress_study.txt')
    Tname = stressname;
end

% Maximum temperature when the value stays the same. Set these to NaN
C(isnan(dz)) = nan;
for i = 1:length(paramv)
    j0(i) = find(dz(:,i)==min(dz(:,i)));
end
plot_tikz(Tv, C, paramv, Tname, Cname, paramname, [fname '.C-T.tex']);
plot_tikz(Tv, dz, paramv, Tname, dzname, paramname, [fname '.dz-T.tex']);
plot_tikz(dz, C, paramv, dzname, Cname, paramname, [fname '.C-dz.tex']);

% Get p-values for linearity of C-T. We should be able to identify when
% the relationship becomes significantly non-linear.
alpha = 0.05;
p_CT = zeros(size(C));
p_dzT = zeros(size(dz));
p_Cdz = zeros(size(C));
for i = 1:length(paramv)
    for j = 1:length(Tv)
        if isnan(C(j,i))
            p_CT = nan;
            p_dzT = nan;
            p_Cdz = nan;
        end
        if j < j0(i)
            x_C = C(j:j0(i),i) - C(j0(i),i);
            x_dz = dz(j:j0(i),i) - dz(j0(i),i);
            x_T = Tv(j:j0(i))' - Tv(j0(i));
        elseif j > j0(i)
            x_C = C(j0(i):j,i) - C(j0(i),i);
            x_dz = dz(j0(i):j,i) - dz(j0(i),i);
            x_T = Tv(j0(i):j)' - Tv(j0(i));
        else
            p(j,i) = 0;
            continue
        end
        % Perform regression to get p-value out
        [~,~,~,~,stats] = regress(x_C, x_T);
        p_CT(j,i) = stats(end);
        [~,~,~,~,stats] = regress(x_dz, x_T);
        p_dzT(j,i) = stats(end);
        [~,~,~,~,stats] = regress(x_C, x_dz);
        p_Cdz(j,i) = stats(end);
    end
end
plot_tikz(Tv, p_CT, paramv, Tname, 'p for C--T', paramname,...
    [fname '.linearity_C-T.tex'], @plot_alpha_cb);
plot_tikz(Tv, p_dzT, paramv, Tname, 'p for dz--T', paramname,...
    [fname '.linearity_dz-T.tex'], @plot_alpha_cb);
plot_tikz(dz, p_Cdz, paramv, dzname, 'p for C--dz', paramname,...
    [fname '.linearity_C-dz.tex'], @plot_alpha_cb);
% Callback to plot alpha level
    function plot_alpha_cb()
        h = ishold();
        hold('on');
        plot(Tv, alpha*ones(size(Tv)), '--k')
        ax = gca();
        if ~ismember(alpha, ax.YTick)
            ax.YTick = sort([ax.YTick alpha]);
        end
        ax.YTickLabel{ax.YTick==0.05} = '(\alpha) 0.05';
        if ~h
            hold('off')
        end
    end

% Compute sensitivity (gradient) of linear region
S_CT = nan(size(paramv));
S_Cdz = nan(size(paramv));
S_dzT = nan(size(paramv));
for i = 1:length(paramv)
    % C-T sensitivity
    % Start of linear region
    j_min = find(p_CT(1:j0(i),i)<alpha, 1, 'first');
    if isempty(j_min)
        j_min = 1;
    end
    % End of linear region
    j_max = find(p_CT(j0(i):end,i)<alpha, 1, 'last');
    if isempty(j_max)
        j_max = length(Tv);
    end
    S_CT(i) = mean(diff(C(j_min:j_max,i))./diff(Tv(j_min:j_max)));
    
    % dz-T sensitivity
    j_min = find(p_dzT(1:j0(i),i)<alpha, 1, 'first');
    if isempty(j_min)
        j_min = 1;
    end
    j_max = find(p_dzT(j0(i):end,i)<alpha, 1, 'last');
    if isempty(j_max)
        j_max = length(Tv);
    end
    S_dzT(i) = mean(diff(dz(j_min:j_max,i))./diff(Tv(j_min:j_max)));
    
    % C-dz sensitivity
    j_min = find(p_Cdz(1:j0(i),i)<alpha, 1, 'first');
    if isempty(j_min)
        j_min = 1;
    end
    % End of linear region
    j_max = find(p_Cdz(j0(i):end,i)<alpha, 1, 'last');
    if isempty(j_max)
        j_max = length(Tv);
    end
    S_Cdz(i) = mean(diff(dz(j_min:j_max,i))./diff(dz(j_min:j_max,i)));
end
plot_tikz(paramv, S_CT, 'C-T sensitivity', {}, paramname, yname, '',...
    [fname 'sensitivity_C-T.tex']);
plot_tikz(paramv, S_dzT, 'dz-T sensitivity', {}, paramname, yname, '',...
    [fname 'sensitivity_dz-T.tex']);
plot_tikz(paramv, S_Cdz, 'C-dz sensitivity', {}, paramname, yname, '',...
    [fname 'sensitivity_C-dz.tex']);

end

function [paramv, Tv, C, dz, cte, paramname] = parse_fname(fname)
%DESCRIPTION
%   param, T, C, dz, cte
%
%INPUTS
%   fname - COMSOL results .txt file to parse

% MATLAB's builtin importdata is pretty effective.
file = importdata(fname);
if iscell(file)
    for i = 1:size(file, 1)
        if file{i,1}(1) ~= '%'
            break
        end
    end
    s = struct();
    s.textdata = file(1:i-1,1);
    s.data = cellfun(@str2double, file(i:end,:));
    file = s;
end

% Get data
data = file.data;

% Parse header.
header = cellfun(@(x)(x(3:end)), file.textdata, 'UniformOutput', false);
headings = strsplit(...
    header{end}, '\s\s+', 'DelimiterType', 'RegularExpression');
header = cellfun(...
    @(x)strsplit(x, ':\s+', 'DelimiterType', 'RegularExpression'),...
    header(1:end-1),...
    'UniformOutput', false);
header = vertcat(header{:})';
header = struct(header{:});
disp(header)
paramname = headings{1};

% Format data two independent vectors and two dependent grids.
Tv = data(:,2)';
for i = 2:size(data,1)
    if ismember(Tv(i), Tv(1:i-1))
        Tv = Tv(1:i-1);
        break
    end
end
paramv = data(1:length(Tv):end,1);
C = reshape(data(:,3), length(Tv), length(paramv));
dz = reshape(data(:,4), length(Tv), length(paramv));
cte = reshape(data(:,5), length(Tv), length(paramv));
end

function plot_tikz(x, y, paramv, xname, yname, paramname, fname, cb)
%DESCRIPTION
%
%
%INPUTS
%
figure()
plot(x, y);
if ~isvector(y)
    if isnumeric(paramv)
        paramv = arrayfun(@(x)num2str(x, [paramname ' %g']), paramv,...
                          'UniformOutput', false);
    end
    legend(paramv);
end
ylabel(yname);
xlabel(xname);
if exist('cb', 'var')
    cb()
end
matlab2tikz(fname);
end