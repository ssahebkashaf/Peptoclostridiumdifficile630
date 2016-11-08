% FLUX_BALANCE Flux-balance analysis of FBA model
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL) performs a basic flux-balance
%    analysis of the FBA model FBAMODEL and returns a biomass-maximizing
%    flux distribution in the vector V.  The maximimum and minimum 
%    synthetic objective possible in a biomass-maximizing flux distribution
%    is given in FMAX and FMIN, respectively.
%
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL, QUIET) performs the analysis
%    and supresses screen output if QUIET is set to true.

function [v, fmax, fmin] = flux_balance(fbamodel, quiet)

if nargin < 2
    quiet = false;
end

param.tmlim  = -1;
param.msglev = 1;
param.save   = 0;

nrxn   = fbamodel.nrxn;
nmetab = fbamodel.nmetab;

yt = fbamodel.present;

A = [ fbamodel.S; 
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn) ];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1) ];
ctype = char('S' * ones(1, nmetab + nnz(~yt)));
vartype = char('C' * ones(1, nrxn));
[v, vbiomass] = glpk(fbamodel.f, A, b, fbamodel.vmin, fbamodel.vmax, ctype, vartype, -1, param);

A = [ fbamodel.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      fbamodel.f' ];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass ];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 1));
[v, vbiomass] = glpk(fbamodel.f, A, b, fbamodel.vmin, fbamodel.vmax, ctype, vartype, -1, param);

[v, vinta] = glpk(fbamodel.g, A, b, fbamodel.vmin, fbamodel.vmax, ctype, vartype, 1, param);

if ~quiet
    fprintf('Biomass flux:    %f\n', fbamodel.f' * v)
    fprintf('Synthetic flux:  [%f, %f]\n', fmin, fmax)
end