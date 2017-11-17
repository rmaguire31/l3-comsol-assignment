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
[results, xv, yv, headings] = parse_fname(fname)
end

function [results, xv, yv, headings] = parse_fname(fname)
%DESCRIPTION
%
%
%INPUTS
%   fname - COMSOL results .txt file to parse

% MATLAB's builtin importdata is pretty effective.
file = importdata(fname)
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
disp(headings)

% Format data two independent vectors and a dependent grid.
yv = nan(1,size(data,1));
for i = 1:size(data,1)
    if ismember(data(i,2), yv)
        yv = yv(1:i-1);
        break
    else
        yv(i) = data(i,2);
    end
end
xv = data(1:length(yv):end,1)';
results = reshape(data(:,end), length(xv), length(yv));
end