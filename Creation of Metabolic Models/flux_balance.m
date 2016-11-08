% FLUX_BALANCE Flux-balance analysis of FBA model
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL) performs a basic flux-balance
%    analysis of the FBA model FBAMODEL and returns a biomass-maximizing
%    flux distribution in the vector V.  The maximimum and minimum 
%    synthetic objective possible in a biomass-maximizing flux distribution
%    is given in FMAX and FMIN, respectively.
%
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL, QUIET) performs the analysis
%    and supresses screen output if QUIET is set to true.

function [v, v1max, v1min, fmax, fmin] = flux_balance(fbamodel, quiet)

if nargin < 2
    quiet = false;
end

param.tmlim  = -1;
param.msglev = 1;
param.save   = 0;

nrxn   = size(fbamodel.rxns,1);
nmetab = size(fbamodel.mets,1);

yt = ones(nrxn,1); %the old yt = fbamodel.present, meaning that all the reactions must be considered as active in the model

A = [ fbamodel.S; 
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn) ];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1) ];
ctype = char('S' * ones(1, nmetab + nnz(~yt)));
vartype = char('C' * ones(1, nrxn));
[v, vbiomass] = glpk(fbamodel.f, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, -1, param); %maximization of fbamodel.f
%[v, vbiomass] = glpk(fbamodel.f, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, 1, param); %minimization of fbamodel.f

A = [ fbamodel.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      fbamodel.f' ];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass ];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 1));
[v1max, fmax] = glpk(fbamodel.g, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, -1, param);
% the "v" that will be considered as final is only the last one (i.e. below this comment). So, for exampla, if we are interested in the flux distribution v associated with the fmin, the [v, fmin] line must be the last one (and it must come after the [v, fmax] line
[v1min, fmin] = glpk(fbamodel.g, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, 1, param);

if ~quiet
    fprintf('Biomass flux:    %f\n', fbamodel.f' * v);
    fprintf('Synthetic flux:  [%f, %f]\n', fmin, fmax);
end