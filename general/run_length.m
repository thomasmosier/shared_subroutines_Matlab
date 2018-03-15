function I = run_length(A)
%Based on "runindex" from Matlab FileExchange by Jos


%Strart of runs:
Q = [false ; diff(A(:)) ~= 0] ;
Rstart = [1 ; find(Q)] ;
I = ones(size(Q)) ;
% cumsum(I) will now increase from 1 to N (= number of elements in V)
% Reset the positions where V changes (Q) by substracting the length of the
% run before it.
I(Q) = 1 - diff(Rstart) ;
% Now we can apply cumsum
I = cumsum(I);