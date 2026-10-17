function M = read_index(f)
%READ_INDEX  index of a KFS monotonic series: one row per test, [c1, p0, e0].
% c1 is the axial strain at the end of the test in TMD and has no meaning in TMU/TMU2 (a copy of the
% TMD column); only p0 and e0 are used. The files use ';' (TMU) or ',' (TMD) as separator.
fid = fopen(f, 'r');
assert(fid > 0, 'read_index: cannot open %s', f);
C = textscan(fid, '%f%f%f', 'Delimiter', {';' ','}, 'CollectOutput', true, 'MultipleDelimsAsOne', false);
fclose(fid);
M = C{1};
M = M(all(isfinite(M), 2), :);
assert(~isempty(M) && size(M, 2) == 3, 'read_index: %s does not have 3 columns.', f);
end
