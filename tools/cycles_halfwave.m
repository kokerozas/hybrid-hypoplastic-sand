function Nvec = cycles_halfwave(q)
%CYCLES_HALFWAVE  cycle axis N of a cyclic record: N advances 0.5 between two reversals of q.
%
% A reversal is a change of sign of q after a stretch whose maximum |q| reaches 30 % of the maximum
% of the record (filters noise around q = 0 in the measured curves). N grows linearly inside each
% half cycle. Used for the KFS numbers of the thesis and for the measured KFS records.
q = q(:); n = numel(q);
thr = 0.3*max(abs(q));
if thr <= 0, Nvec = zeros(n,1); return; end
s = sign(q); s(s==0) = 1;
raw = find(diff(s) ~= 0);
keep = zeros(numel(raw),1); nk = 0; runmax = 0; prev = 1;
for j = 1:numel(raw)
    idx = raw(j);
    runmax = max(runmax, max(abs(q(prev:idx))));
    prev = idx + 1;
    if runmax >= thr, nk = nk + 1; keep(nk) = idx; runmax = 0; end
end
cr = [1; keep(1:nk); n];
Nvec = zeros(n,1); cval = 0;
for i = 1:numel(cr)-1
    seg = cr(i):cr(i+1);
    Nvec(seg) = linspace(cval, cval+0.5, numel(seg));
    cval = cval + 0.5;
end
end
