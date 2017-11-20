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
fname = fname(1:end-4);

Tname = 'Temperature \si{\celsius}';
Cname = 'Capacitance \si{\femto\farad}';
dzname = 'Displacement \si{\micro\meter}';
if strcmp(fname, 'material_study')
    paramname = 'Thermal expansion coefficient \si\{\per\kelvin}';
    material = {'\ce{Au}', '\ce{Al}', '\ce{Ti}'};
    paramv = arrayfun(...
        @(i,x)[material{i} num2str(x,'\\tiny($\\alpha=\\SI{%g}{\\per\\kelvin}$)')], ...
        paramv, cte(1,:)', 'UniformOutput', false);
end

% Maximum temperature when the value stays the same. Set these to NaN
C(isnan(dz)) = nan;
plot_tikz(Tv, C, paramv, Tname, Cname, paramname, [fname '.C-T.tex']);
plot_tikz(Tv, dz, paramv, Tname, dzname, paramname, [fname '.dz-T.tex']);
plot_tikz(dz, C, paramv, dzname, Cname, paramname, [fname '.C-dz.tex']);

% Get p-values for linearity of C-T. We should be able to identify when
% the relationship becomes significantly non-linear.
    function [p, F, b] = regression(X,Y)
        p = zeros(length(X),length(X),size(Y,2));
        F = zeros(length(X),length(X),size(Y,2));
        b = cell(length(X),length(X),size(Y,2));
        for ii = 1:size(Y,2)
            for jj = 1:length(X)
                for kk = jj:length(X)
                    if any(isnan(Y(jj:kk,ii)))
                        p([kk,jj],[jj,kk],ii) = nan;
                        continue
                    end
                    if isvector(X)
                        x = reshape(X(jj:kk), [], 1);
                    else
                        x = X(jj:kk,ii);
                    end
                    % Perform regression to get p-values out
                    [beta,~,~,stats] = regress(...
                          Y(jj:kk,ii), [ones(1+kk-jj,1) x]);
                    b([kk,jj],[jj,kk],ii) = {beta};
                    F([kk,jj],[jj,kk],ii) = stats(end-1);
                    p([kk,jj],[jj,kk],ii) = stats(end);
                end
            end
        end
    end
[p_CT, F_CT, b_CT] = regression(Tv,C);
[p_dzT, F_dzT, b_dzT] = regression(Tv,dz);
[p_Cdz, F_Cdz, b_Cdz] = regression(dz(:),C(:));

disp('C-T')
disp(p_CT)
disp('dz-T')
disp(p_dzT)
disp('C-dz')
disp(p_Cdz)

% Find linear region
    function lin = linear_region(p)
        alpha = 0.05;
        lin = cell(size(p,3),3);
        for ii = 1:size(p,3)
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
lin_Cdz = linear_region(p_Cdz);

disp('C-T')
disp(lin_CT)
disp('dz-T')
disp(lin_dzT)
disp('C-dz')
disp(lin_Cdz)

% Calculate the sensitivity in the linear region using our regression
% matrix
    function S = sensitivity(lin, b)
        S = ones(size(lin,1),1);
        for ii = 1:size(lin,1)
            beta = b{lin{ii,:}};
            S(ii) = beta(end);
        end
    end
S_CT = sensitivity(lin_CT, b_CT);
S_dzT = sensitivity(lin_dzT, b_dzT);
S_Cdz = sensitivity(lin_Cdz, b_Cdz);

disp('C-T')
disp(S_CT)
disp('dz-T')
disp(S_dzT)
disp('C-dz')
disp(S_Cdz)


% Get points for the best fit line in the linear region using our 
% regression matrix
    function [x, y] = bestfit(lin, b, X)
        x = ones(2,size(lin,1));
        y = ones(2,size(lin,1));
        for ii = 1:size(lin,1)
            if isvector(X)
                [kk,jj] = lin{ii,1:2};
                x(:,ii) = X([jj kk]);
            else
                [kk,jj] = lin{ii,1:2};
                x(:,ii) = X([jj kk], ii);
            end
            beta = b{lin{ii,:}};
            y(:,ii) = [ones(2,1) x(:,ii)]*beta;
        end
    end
[x_CT,y_CT] = bestfit(lin_CT, b_CT, Tv);
[x_dzT,y_dzT] = bestfit(lin_dzT, b_dzT, Tv);
[x_Cdz,y_Cdz] = bestfit(lin_Cdz, b_Cdz, dz(:));

plot_tikz(Tv, C, paramv, Tname, Cname, paramname,...
     [fname '.linearity_C-T.tex'], @()plot_bestfit_cb(x_CT,y_CT));
plot_tikz(Tv, dz, paramv, Tname, dzname, paramname,...
    [fname '.linearity_dz-T.tex'], @()plot_bestfit_cb(x_dzT,y_dzT));
[~,idx] = sort(dz(:));
plot_tikz(dz(idx), C(idx), {}, dzname, Cname, '',...
    [fname '.linearity_C-dz.tex'], @()plot_bestfit_cb(x_Cdz,y_Cdz));
    function plot_bestfit_cb(x,y)
        h = ishold();
        hold('on');
        for ii = 1:size(x,2)
            p = plot(x(:,ii), y(:,ii), '--');
            p.LineWidth = 4*p.LineWidth;
        end
        if ~h
            hold('off')
        end
    end
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
% Two-fold symmetry => 4*C
C = 4*reshape(data(:,3), length(Tv), length(paramv));
if (size(data,2) >= 4)
    dz = reshape(data(:,4), length(Tv), length(paramv));
else
    dz = C;
end
if (size(data,2) >= 4)
    cte = reshape(data(:,5), length(Tv), length(paramv));
else
    cte = dz;
end
end

function plot_tikz(x, y, paramv, xname, yname, paramname, fname, cb)
%DESCRIPTION
%
%
%INPUTS
%
markers = {
    '-+','-o','-^','-s','-p','-h','-v','-d'};
figure()
hold('on');
for i = 1:size(y,2)
    if isvector(x)
        plot(x, y(:,i), markers{i});
    else
        plot(x(:,i), y(:,i), markers{i});
    end
end
hold('off');

if ~isvector(y)
    if isnumeric(paramv)
        paramv = arrayfun(@(x)[paramname num2str(x,' %g')], paramv,...
                          'UniformOutput', false);
    end
   legend(paramv, 'AutoUpdate', 'off', 'Location', 'northoutside');
end
ylabel(yname);
xlabel(xname);

% Callback if we need to add anything else.
if exist('cb', 'var')
    cb()
end
matlab2tikz(['tex/tikz/' fname], ...
            'parseStrings', false, ...
            'width', '0.4\textwidth');
end