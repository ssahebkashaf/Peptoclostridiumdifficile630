

warning('off', 'all')

num_cores = 4;
%addAttachedFiles(gcp,{'glpk.m','glpkcc'});
expFBA(140,6000,'iCdiff834',4995,num_cores);

