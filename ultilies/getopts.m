%% ------------------SUBFUNCTION-----------------------------
function opts = getopts(opts)

if ~isfield(opts,'maxitr')
    opts.maxitr = 5;
end
if ~isfield(opts,'beta1')
    opts.beta1 = 5;
end
if ~isfield(opts,'beta2')
    opts.beta2 = 20;
end

if ~isfield(opts,'gamma')
    opts.gamma = 1.618;
end
if ~isfield(opts,'relchg')
    opts.relchg = 1.e-3;
end
if ~isfield(opts,'print')
    opts.print = 0;
end
