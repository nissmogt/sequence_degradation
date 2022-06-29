% Copyright 2014 - by Magnus Ekeberg (magnus.ekeberg@gmail.com)
% All rights reserved
% 
% Permission is granted for anyone to copy, use, or modify this
% software for any uncommercial purposes, provided this copyright 
% notice is retained, and note is made of any changes that have 
% been made. This software is distributed without any warranty, 
% express or implied. In no event shall the author or contributors be 
% liable for any damage arising out of the use of this software.
% 
% The publication of research using this software, modified or not, must include 
% appropriate citations to:
%
% 	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
% 	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013) 
%
%	M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
%	maximization for direct-coupling analysis of protein structure
%	from many homologous amino-acid sequences, arXiv:1401.4832
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [fval,grad] = g_r(wr,Y,weights,N,q,lambdah,halflambdaJ,r)
%Evaluates (regularized) g_r using the mex-file.
h_r=reshape(wr(1:q),1,q);
J_r=reshape(wr(q+1:end),q,q,N-1);

r=int32(r);
[fval,grad1,grad2] = g_rC(Y-1,weights,h_r,J_r,[lambdah;halflambdaJ],r);
grad = [grad1(:);grad2(:)];


