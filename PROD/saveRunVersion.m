function runDir = saveRunVersion(scriptName, params, outputFiles, label)
%SAVERUNVERSION Αποθηκεύει τα αποτελέσματα ενός run σε versioned φάκελο μαζί
% με πλήρες snapshot παραμέτρων και diff vs την προηγούμενη εκτέλεση.
%
%   runDir = saveRunVersion(scriptName, params, outputFiles, label)
%
% scriptName  : char, π.χ. 'test_simulation'
% params      : struct με τις παραμέτρους του run που θέλουμε να καταγραφούν
% outputFiles : cellstr με πλήρεις διαδρομές αρχείων που παρήγαγε το run
%               (αντιγράφονται στον φάκελο)· {} ή [] αν κανένα
% label       : char, προαιρετικό tag στο όνομα φακέλου· '' αν κανένα
%
% Δομή: Results/runs/<YYYYMMDD_HHMMSS>_<script>[_<label>]/
%   params.mat   ακριβές restore
%   params.txt   αναγνώσιμο snapshot (key = value)
%   changed.txt  diff παραμέτρων vs την προηγούμενη εκτέλεση του ίδιου script
%   <CSV/PNG του run>

if nargin < 3 || isempty(outputFiles), outputFiles = {}; end
if nargin < 4, label = ''; end
if ischar(outputFiles), outputFiles = {outputFiles}; end

rootDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results', 'runs');
if ~isfolder(rootDir), mkdir(rootDir); end

stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
name  = [stamp '_' scriptName];
if ~isempty(label), name = [name '_' label]; end
runDir = fullfile(rootDir, name);
mkdir(runDir);

% -- Μεταδεδομένα περιβάλλοντος --
params.meta = struct('timestamp', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ... %#ok<TNOW1,DATST>
                     'script', scriptName, 'label', label, ...
                     'matlabVersion', version, 'gitCommit', localGitCommit());

% -- Snapshot --
save(fullfile(runDir, 'params.mat'), 'params');
paramLines = flattenStruct(params, '');
writeLines(fullfile(runDir, 'params.txt'), paramLines);

% -- Diff vs προηγούμενη εκτέλεση του ίδιου script --
prev = findPreviousRun(rootDir, scriptName, runDir);
if isempty(prev)
    writeLines(fullfile(runDir, 'changed.txt'), ...
        {sprintf('Καμία προηγούμενη εκτέλεση του %s - baseline run.', scriptName)});
else
    d = diffKeyValue(readLines(fullfile(prev, 'params.txt')), paramLines);
    [~, prevName] = fileparts(prev);
    header = {sprintf('Σύγκριση vs %s', prevName), ''};
    if isempty(d), d = {'(καμία αλλαγή παραμέτρων)'}; end
    writeLines(fullfile(runDir, 'changed.txt'), [header, d]);
end

% -- Αντιγραφή outputs --
for i = 1:numel(outputFiles)
    if ~isempty(outputFiles{i}) && isfile(outputFiles{i})
        copyfile(outputFiles{i}, runDir);
    end
end

fprintf('Run version -> %s\n', runDir);
end

% ======================================================================

function c = localGitCommit()
c = 'unknown';
try
    here = fileparts(mfilename('fullpath'));
    [st, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', here));
    if st == 0
        c = strtrim(out);
        [st2, out2] = system(sprintf('git -C "%s" status --porcelain', here));
        if st2 == 0 && ~isempty(strtrim(out2))
            c = [c ' (dirty)'];
        end
    end
catch
end
end

function lines = flattenStruct(s, prefix)
lines = {};
if ~isstruct(s), return; end
f = fieldnames(s);
for i = 1:numel(f)
    key = f{i};
    val = s.(key);
    if isempty(prefix), fullkey = key; else, fullkey = [prefix '.' key]; end
    if isstruct(val)
        lines = [lines, flattenStruct(val, fullkey)]; %#ok<AGROW>
    elseif isnumeric(val) || islogical(val)
        lines{end+1} = sprintf('%s = %s', fullkey, mat2str(val, 8)); %#ok<AGROW>
    elseif ischar(val)
        lines{end+1} = sprintf('%s = %s', fullkey, val); %#ok<AGROW>
    elseif isstring(val)
        lines{end+1} = sprintf('%s = %s', fullkey, strjoin(cellstr(val(:).'), ', ')); %#ok<AGROW>
    else
        lines{end+1} = sprintf('%s = <%s>', fullkey, class(val)); %#ok<AGROW>
    end
end
end

function prev = findPreviousRun(rootDir, scriptName, currentRunDir)
prev = '';
d = dir(fullfile(rootDir, ['*_' scriptName '*']));
d = d([d.isdir]);
names = sort({d.name});
[~, cur] = fileparts(currentRunDir);
names = names(~strcmp(names, cur));
if ~isempty(names)
    prev = fullfile(rootDir, names{end});
end
end

function writeLines(path, lines)
fid = fopen(path, 'w', 'n', 'UTF-8');
if fid == -1, error('saveRunVersion:write', 'Αδυναμία εγγραφής %s', path); end
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
fclose(fid);
end

function lines = readLines(path)
lines = {};
if ~isfile(path), return; end
fid = fopen(path, 'r', 'n', 'UTF-8');
raw = fread(fid, '*char').';
fclose(fid);
lines = strsplit(raw, newline);
lines = lines(~cellfun(@isempty, lines));
end

function d = diffKeyValue(oldLines, newLines)
oldMap = kvMap(oldLines);
newMap = kvMap(newLines);
d = {};
keys = union(oldMap.keys, newMap.keys);
for i = 1:numel(keys)
    k = keys{i};
    if startsWith(k, 'meta.'), continue; end   % timestamp/commit αλλάζουν πάντα
    hasOld = isKey(oldMap, k);
    hasNew = isKey(newMap, k);
    if hasOld && hasNew
        if ~strcmp(oldMap(k), newMap(k))
            d{end+1} = sprintf('  ~ %s : %s  ->  %s', k, oldMap(k), newMap(k)); %#ok<AGROW>
        end
    elseif hasNew
        d{end+1} = sprintf('  + %s = %s', k, newMap(k)); %#ok<AGROW>
    else
        d{end+1} = sprintf('  - %s  (ήταν: %s)', k, oldMap(k)); %#ok<AGROW>
    end
end
end

function m = kvMap(lines)
m = containers.Map('KeyType', 'char', 'ValueType', 'char');
for i = 1:numel(lines)
    t = regexp(lines{i}, '^(.*?)\s*=\s*(.*)$', 'tokens', 'once');
    if numel(t) == 2
        m(strtrim(t{1})) = strtrim(t{2});
    end
end
end
