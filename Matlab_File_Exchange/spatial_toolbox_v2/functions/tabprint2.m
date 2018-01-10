function tabout=tabprint2(z, info)
%
% PURPOSE: puts an (n_rows x n_cols) matrix in string cell array tabular format with decimal alignment
%                     and commas every three digits in math mode. One can print this by 
%                     simply listing the name of the cell array in a program or in the command line (e.g., tabout).
%                     Using format long will usually lead Matlab to print all the column headings for cell arrays.
%                     Optionally allows control of precision. 
%                     Optionally allows for column and row headings.  
%                     Optionally places table title.  When title is present, the routine directly prints the table with title.
%                     In the case with the title, one may wish to suppress the printing by the structure name by giving a semicolon
%                     after the function invocation.
%                     
%---------------------------------------------------
% USAGE:    tabprint2(z, info)
% where:
%        z         = (n_rows x n_cols) matrix (or vector) to be printed. Assume this contains real numbers.
%        info      = a structure containing printing options
%        info.cnames = an (n_cols x 1) string vector of names for columns (optional)
%                      e.g. info.cnames = strvcat('col1','col2');
%                      (default = integer column headings)
%        info.rnames = an (n_rows+1 x 1) string vector of names for rows (optional)
%                      e.g. info.rnames = strvcat('Rows','row1','row2');
%                      (default = integer row labels)
%        info.fmt    = a format string, e.g., '%12.6f' or '%12d' (default = '%25.4f'). Maximum of 9 decimal places. 
%        info.caption = a string giving the caption or title of the table (default is a blank)
%                    
%---------------------------------------------------
% EXAMPLES:
%        in.cnames = strvcat('col1','col2');
%        in.rnames = strvcat('rowlabel','row1','row2');
%        in.caption=['Estimates'];
%   
%        tabprint2(y, in), prints entire matrix, column headings, row labels, and title of the table
%
%  or: supplying in.fmt with just  one the column, row label fields will result
%        in integer labels for the missing field(s).
%
%  or: tabprint2(y) prints entire matrix with integer column headings and integer row labels 
%
%  Try tabprint2(randn(10,3)*100000) to see how it works
%
%NOTES:
%
% tabprint2 is partially based upon mprint and lprint written by James P. LeSage, www.spatial-econometrics.com
% and on int2str2 written by Olivier Ledoit. It uses the same syntax as mprint and lprint with the addition of caption and label
% fields. It is compatible with and derived from the function latextab2 (latex table printing).
%
%Written by Kelley Pace, 1/1/03, www.spatial-statistics.com
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Input Checking and Default Setting Section %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%checks to see if there are just one or two arguments
error(nargchk(1,2,nargin))

%checking to see if the first argument is numeric
if isnumeric(z)==0
    error('z must be a numeric matrix')
end

%creating a default structure when none exists
if nargin==1
      info.default=1;%have to create some field to make structure non-empty     
end

%testing for the presence of a structure  when second argument is present  
if nargin==2
  if ~isstruct(info)
    error('If you supply options, you must supply the options as a structure variable'); 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sets the format, when not supplied

if isfield(info,'fmt')==0
    info.fmt='%25.4f';
end

if isstr(info.fmt)
%this assumes that there is 9 or fewer digits of precision
dposition=regexpi(info.fmt,'\.')+1;
dprecision=str2num(info.fmt(dposition));
else
error('info.fmt should be a string of form %integer.integerf   (eg.  %15.4f) -- also require fewer than 9 decimal places.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%need some dimensions

[n_rows, n_cols]=size(z);
n_rows_table=n_rows+1;
n_cols_table=n_cols+1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sets row names 

if isfield(info,'rnames')==1
v1=cellstr(info.rnames); 

%checks to see if the supplied row names matches the number required
if length(v1)~=n_rows_table
    error('error tprint: The number of row labels is incorrect')
end;

else
%sets up some default row labels
for ijij=1:n_rows_table
    v1(ijij)=cellstr(num2str(ijij-1));
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sets column headings

 if isfield(info,'cnames')==1
v2=cellstr(info.cnames);   

%checks to see if the supplied column names matches the number required
if length(v2)~=n_cols
    error('error tprint: The number of column labels is incorrect')
end;

else
%sets up some default column headings
for ijij=1:n_cols
    v2(ijij)=cellstr(num2str(ijij));
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Setting up cells with the right format %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%turns all contents into a vector
zall=z(:);

%integer part of number
zall_left=fix(zall); 

%for each number, convert it into a string under control of info.fmt. Take the integer part of number and 
%separately convert it to a string with commas every three numbers. Merge the decimal part with this.
%Store the combined string in the appropriate cell.
for iii=1:length(zall)
    stemp1=num2str(zall(iii), info.fmt);
    stemp2=int2str2(zall_left(iii));
    stemp3=[stemp2 '.' stemp1((end-dprecision+1):end) ];
    zs1(iii)=cellstr(stemp3);
end;

%reshape the numerical part back
zs2=reshape(zs1, n_rows, n_cols);

%preallocating the output cells
vc=cell(n_rows_table, n_cols_table);

%putting together the row, column headings, and the numerical contents
vc(:,1)=v1;
vc(1,2:end)=v2';
vc(2:end,2:end)=zs2;

%converting back to strings, right justifying, and converting back to cell
vc2=cellstr(strjust(char(vc),'right'));

%reshaping back to table dimensions
tabout=reshape(vc2, n_rows_table, n_cols_table);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%printing directly to screen when caption is present

if isfield(info,'caption')==1
    %sets off by a couple of blank lines
    disp(blanks(2)');
    %display title left justified
    disp(info.caption);
    %set off title slightly from variable name
    disp(' ');
    %reduce line spacing between variable name and table
    format compact
    %give long format to display more cell contents (for long column headings)
     format long;
     %invoke the matlab cell print routine
    tabout
    %set off by a final blank line
    disp(' ');
    %restore to a default format
    format short;
    format loose;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Integer to strings with commas subfuction from Olivier Ledoit %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function s=int2str2(n)
%
% converts integer n into string s
% inserts commas every three digits
% n must be a scalar
%
% Example: disp(int2str2(-1e6))
% 
% Author: Olivier Ledoit
%         olivier.ledoit@csfb.com

function s=int2str2(n)

s=int2str(abs(n));
i=mod(-length(s),3);
s=[repmat('0',[1 i]) s];
j=length(s)/3;
s=reshape([reshape(s,[3 j]);repmat(',',[1 j])],[1 4*j]);
s([1:i 4*j])=[];
s=[repmat('-',[1 n<0]) s];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%