% ****************************************************************
% kernpdf_modify.m
%
% This function creates a kernel-based empirical probability distribution
%   function using a normal kernel. It is passed a 1-D array of data
%   points, along with an array of quantiles to calculate,
%   and returns the values at the specified quantiles
%
% Note: data pt "0" is treated differently -- using "delta" function here
%   i.e. histogram approach. Also, kernels are rescaled if their probability
%   significantly leaks over into the negative axis.
%
% The bandwith of the kernel is (was?) prescribed from the reference:
%   Silverman, B. W. (1986), Density Estimation for Statistics and Data
%   Analysis, New York: Chapman and Hall.
%
% Assuming past quantile values (quans) are in ascending order (sorted)
% *****************************************************************
%  HRC on 11/7/2016

function quanval = kernpdf_modify(data,quans)

KNORM=1 ;

SKEWDAT=1 ;
if SKEWDAT == 1
    PWR=1.8 ; % fractional value to change xdel
end

% plot data and generate historgram, PLTDAT=1;
PLTDAT=0 ;


% constants

n = size(data,1) ; % number of array elements
XMAX=1.3*max(data(:))  ; % largest value to evaluate PDF at
XMIN = 0. ; % minimum value to evaluate PDF at
NQUAN= length(quans(:)) ; % number of quantiles
NEVAL=2*NQUAN ; % number of evaluation points
MINPPT=1.e-5  ; % minimum precip value to group as 'no precip' [m/day]
RDEL=1/2.^32 ;  % extra rounding error "slop" added to bounds

pdfarr=zeros(NEVAL,1) ; % calculated empirical PDF

% first, sort data

sdata=sort(data(:)) ;
data = sdata ;

% next, determine start of data greater than min threshold

idx = find(data >  MINPPT) ;
if length(idx) < NQUAN
    quanval = zeros(NQUAN,1) ;
    return
end
if length(idx) == 0
    quanval = zeros(NQUAN,1) ;
%   sprintf('no data above threshold')
    return ;
end
inon0 = idx(1) ;
% sprintf('first data pt above threshold at: %d',inon0) 

delx=zeros(NEVAL,1) ;
if SKEWDAT == 0
    delx(1)=2.*MINPPT ;
    delx(2:NEVAL)=(XMAX-MINPPT)/(NEVAL-1.5) ;
end
if SKEWDAT == 1
    sknss=skewness(data,1) ;
%   sprintf('sknss=%f',sknss) 
    afac=PWR^(sknss/(NEVAL-1)) ;
%   sprintf('afac = %f',afac)  
    % note: not including first eval pt ("1") in calc
    atot=0.5 ;
    for i=2:NEVAL-2
        atot=atot+afac^i ;
    end
    atot=atot+0.5*afac^(NEVAL-1) ;
%   sprintf('atot=%f',atot) 
    dx0=(XMAX-XMIN-MINPPT)/atot ;
    delx(1) = 2.*MINPPT ;
    for i= 2:NEVAL
        delx(i)=(afac^(i-1))*dx0 ;
    end
end
% sprintf('delx(NEVAL)=%f',delx(NEVAL)) 

% next, evaluation points

xeval=zeros(NEVAL,1) ;
xeval(1)=XMIN ;
for i=2:NEVAL
    xeval(i)=xeval(i-1)+0.5*(delx(i-1)+delx(i)) ;
end
if xeval(NEVAL) ~= XMAX
%   sprintf('xeval(NEVAL) = %f XMAX = %f XMAX-xeval(NEVAL) = %f',xeval(NEVAL),XMAX,XMAX-xeval(NEVAL)) 
end

% lower and upper interval bounds
lbnd=zeros(NEVAL,1) ;
ubnd=zeros(NEVAL,1) ;
lbnd(1) = XMIN ;
ubnd(1) = MINPPT ;  %  note: not used -- included for completeness
for i=2:NEVAL-1
    lbnd(i)=ubnd(i-1) ;
    ubnd(i)=lbnd(i)+delx(i) ;
end
lbnd(NEVAL) = ubnd(NEVAL-1) ;
ubnd(NEVAL) = XMAX ;

RNG=50 ; % data pts (+/- interval) around data pt to calc bandwidth
TOTRNG=2*RNG+1 ;


RESCL=2.3 ;

hband=zeros(n,1) ;
for i=inon0+1:n
    istrt=i-RNG ;
    istop=i+RNG ;
% if istrt LT 1 then istrt=1
    if istrt < inon0 + 1
        istrt=inon0 + 1;
        istop=inon0 + 1 +2*RNG ;
    end
    if istop  > n 
        istop=n ;
        istrt=istop-2*RNG ;
    end
    sig=std(data(istrt:istop));  % standard deviation of the data
    R=iqr(data(istrt:istop));  % inter-quartile range of the data
    dumarr=[sig,R/1.34] ;
%   sprintf('i=%d istrt=%d istop= %d',i,istrt,istop')
%   sprintf('sig= %f R = %f dumarr = %f %f' , sig, R, dumarr(1),dumarr(2)) 
    hband(i)=min(dumarr(:)) ; % bandwidth
%   sprintf('hband(%d)= %f',i,hband(i))
end

hbeval=zeros(n,1) ;
ilow=1 ;
for idat=inon0+1:n
    tttt = true ;
    while tttt
        ilow=ilow+1 ;
        tttt = ~ (data(idat) >= lbnd(ilow) & data(idat) <= ubnd(ilow) ) ;
    end
    hbeval(idat)=delx(ilow) ;
    ilow=ilow-1 ;

    if hbeval(idat) > hband(idat)
        hband(idat)=hbeval(idat) ;
    end
    if data(idat)/hband(idat) < RESCL
        hband(idat)=data(idat)/RESCL ;
%       sprintf('rescaling datapt %d',idat)
    end
end

normc=zeros(n,1) ;
for idat=inon0+1:n
    if KNORM == 0
        normc(idat)=sqrt(2.)/(hband(idat)*sqrt(pi)* ...
          (1.+erf(data(idat)/(hband(idat)*sqrt(2.))))) ;
    end
    if KNORM == 1
        normc(idat)=1./(hband(idat)*sqrt(2.*pi)) ;
    end
end

pdfarr(1)=(inon0+1)/double(n)/(delx(1)/2.) ;

for ieval=2:NEVAL
    xsum=0. ;
% note: only evaluating data GT MINPPT
    for idat=inon0+1:n
        xt=(xeval(ieval) - data(idat))/hband(idat) ;
        if KNORM == 0
            xsum=xsum+normc(idat)*exp((-1/2.)*xt^2) ;
        end
        if KNORM == 1
            xtimg=((xeval(ieval)-MINPPT) + (data(idat)-MINPPT))/hband(idat) ;
            xsum=xsum+normc(idat)*(exp((-1/2.)*xt^2) + exp((-1/2.)*xtimg^2)) ;
        end
    end
% normalize result and assign to PDF array
    pdfarr(ieval)=xsum/double(n) ;
end

cdfarr=pdfarr ; % cumulative probability distribution
cdfarr(1)=pdfarr(1)*delx(1)/2. ;
for ieval=2:NEVAL-1
    cdfarr(ieval)=cdfarr(ieval-1)+pdfarr(ieval)*delx(ieval) ;
end
cdfarr(NEVAL) = cdfarr(NEVAL-1)+pdfarr(NEVAL)*delx(NEVAL)/2. ; % last element
% sprintf('normalization test: %f',cdfarr(NEVAL)) 

% adjust in case not normalized (excluding 1-pt)
normfac=(1.-cdfarr(1))/(cdfarr(NEVAL)-cdfarr(1)) ;
pdfarr(2:NEVAL) = normfac*pdfarr(2:NEVAL) ;
% and now recalc cdfarr w/ updated pdfarr
pdfplt=pdfarr ;
pdfplt(1)=pdfarr(1)*delx(1)/2. ;
pptavg=pdfplt(1)*xeval(1); 
cdfarr(1)=pdfarr(1)*delx(1)/2. ;  % this calc redundant, done for clarity
for ieval=2:NEVAL-1
    pdfplt(ieval)=pdfarr(ieval)*delx(ieval) ;
    pptavg=pptavg+pdfplt(ieval)*xeval(ieval) ;
    cdfarr(ieval)=cdfarr(ieval-1)+pdfarr(ieval)*delx(ieval) ;
end
pdfplt(NEVAL) = pdfarr(NEVAL)*delx(NEVAL)/2. ; % last element
pptavg=pptavg+pdfplt(NEVAL)*xeval(NEVAL) ;
pdfplt = n.*pdfplt ;
cdfarr(NEVAL) = 1. ;

quanval=zeros(NQUAN,1) ; 

ieval=0 ;
for iquan=1:NQUAN
    if quans(iquan) <= cdfarr(1)
        quanval(iquan) = XMIN ;
    else
        tttt = true ;
        while tttt
           ieval=ieval+1 ;
           tttt = ~ (quans(iquan) > cdfarr(ieval) & quans(iquan) <= cdfarr(ieval+1)) ;
        end
%       display(ieval)
% linearly interpolating
        ratio = (quans(iquan)-cdfarr(ieval))/(cdfarr(ieval+1)-cdfarr(ieval)) ;
        quanval(iquan) = xeval(ieval) + ratio*(xeval(ieval+1)-xeval(ieval)) ;
        ieval=ieval-1 ;
    end
end
