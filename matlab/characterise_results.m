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
    j0(i) = find(abs(dz(:,i))==min(abs(dz(:,i))));
end
plot_tikz(Tv, C, paramv, Tname, Cname, paramname, [fname '.C-T.tex']);
plot_tikz(Tv, dz, paramv, Tname, dzname, paramname, [fname '.dz-T.tex']);
plot_tikz(dz, C, paramv, dzname, Cname, paramname, [fname '.C-dz.tex']);

% Get p-values for linearity of C-T. We should be able to identify when
% the relationship becomes significantly non-linear.
    function [p, F, b] = regression(x,y)
        p = zeros(length(x),length(x),length(paramv));
        F = zeros(length(x),length(x),length(paramv));
        b = cell(length(x),length(x),length(paramv));
        for ii = 1:length(paramv)
            for jj = 1:length(x)
                for kk = jj:length(x)
                    if any(isnan(y(jj:kk,ii)))
                        p([kk,jj],[jj,kk],ii) = nan;
                        continue
                    end
                    % Perform regression to get p-values out
                    [beta,~,~,stats] = regress(...
                          x(jj:kk)', [ones(1+kk-jj,1) y(jj:kk,ii)]);
                    b([kk,jj],[jj,kk],ii) = {beta};
                    F([kk,jj],[jj,kk],ii) = stats(end-1);
                    p([kk,jj],[jj,kk],ii) = stats(end);
                end
            end
        end
    end
[p_CT, F_CT, b_CT] = regression(Tv,C);
[p_dzT, F_dzT, b_dzT] = regression(Tv,dz);

disp('C-T')
disp(p_CT)
disp('dz-T')
disp(p_dzT)

% Find linear region
    function lin = linear_region(p)
        alpha = 0.05;
        lin = cell(length(paramv),2);
        for ii = 1:length(paramv)
            % Default linear region
            lin(ii,:) = {1, 1, ii};
            brake = false;
            for diag = size(p,1)-1:-1:0
                for jj = 1:size(p,1)-diag
                    if p(jj,jj+diag,ii) < alpha
                        lin(ii,:) = {jj, jj+diag, ii};
                        brake = true;
                        break
                    end
                end
                if brake
                    break
                end
            end
        end
    end
lin_CT = linear_region(p_CT);
lin_dzT = linear_region(p_dzT);

disp('C-T')
disp(lin_CT)
disp('dz-T')
disp(lin_dzT)

% Calculate the sensitivity in the linear region using regression
    function S = sensitivity(lin, b)
        S = ones(length(paramv),1);
        for ii = 1:length(paramv)
            S(ii) = b{lin{ii,:}};
        end
    end
S_CT = sensitivity(lin_CT, b_CT);
S_dzT = sensitivity(lin_dzT, b_dzT);

disp('C-T')
disp(S_CT)
disp('dz-T')
disp(S_dzT)

return
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
for i = 1:length(paramv)
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