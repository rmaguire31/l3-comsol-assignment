function analyse_prestress(fname)
%CHARACTERISE_RESULTS Plots characteristics of COMSOL results `fname`
%
%DESCRIPTION
%
%
%INPUTS
%   fname - COMSOL results .txt file to parse
%
%COPYRIGHT (C) Russell Maguire 2017
[S0v, Cv, dzv] = parse_fname(fname);
fname = fname(1:end-4);
S0name = 'Pre-stress \si{\giga\pascal}';
dzname = 'Displacement \si{\micro\meter}';

plot(S0v, dzv, '-');
xlabel(S0name)
ylabel(dzname)

hold('on')
plot([S0v(1) S0v(end)], [100 100], '--k');

matlab2tikz(['tex/tikz/' fname], ...
            'parseStrings', false, ...
            'width', '0.35\textwidth');
end

function [tv, Tv, dzv] = parse_fname(fname)
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

% Format data vectors
tv = data(:,1)';
Tv = data(:,2)';
dzv = data(:,3)';
end