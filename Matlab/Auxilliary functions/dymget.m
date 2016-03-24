function data = dymget(dymstr, name)
% dymget  get dymola simulation result for a specific variable 
% 
%    d = dymget(dymstr, name);        dymstr is a structure obtained by executing dymload
%                                     name is the FULL variable name
%                                     d is a vector containing the data
%
%   For simple variables a column vector is returned.
%   For arrays a cell array is returned.
%   For records or sub-models a struct is returned,
%   where each element has the name of the correspond 
%   component/variable.
%
% Issues:
%   If a struct contains a member whose name start
%   with _... the name is changed to FIX_... 
%   since Matlab structs cannot have names starting with _
%
% See also: dymbrowse, dymload

% Originally written by Matthijs Langelaar, DLR.
% Change log:
% 2001-03-09 Hans Olsson, Dynasim AB
%   Restructured code
%   Correct string matching for partially matched names
%   Handles multiple matches correctly
%   Code to handle arrays and structs
%   Performance optimization for string matching
%
% Copyright (c) 2000 DLR,
% Copyright (c) 2001 Dynasim AB, Sweden
%    All rights reserved.

% find index
[index,typ]=myStrMatch(name,dymstr.name);
if length(index)==0 
   error('Invalid variable name');
end
data=dymgetInner(name,index,typ,dymstr,dymstr.name(index,:));

function data=dymgetInner(name,index,typ,dymstr,subnames)
% Handles the actual data
% Recursive
switch(typ)
case 0,
   data=unpackElement(dymstr,index);
case 1,
   % Array: More difficult due to array of components.
   i=1;
   while i<=length(index)
      s=subnames(i,length(name)+2:end);
      I=find(s==']');
      s=s(1:I(1)-1);
      newname=[name,'[',s,']'];
      [index2i,typ2]=myStrMatch(newname,subnames);
      index2=index(index2i);
      value=dymgetInner([name,'[',s,']'],index2,typ2,dymstr,subnames(index2i,:));
      eval(['data{',s,'}=value;']);
      i=i+length(index2);
   end   
case 2,
   % A struct
   i=1;
   first=1;
   while i<=length(index)
      s=subnames(i,length(name)+2:end);
      I=find(((s=='.')|(s=='[')) & cumsum((s=='(')-(s==')'))==0);
      if ~isempty(I)
         s=s(1:I(1)-1);
      end
      newname=[name,'.',s];
      [index2i,typ2]=myStrMatch(newname,subnames);
      index2=index(index2i);
      val=dymgetInner(newname,index2,typ2,dymstr,subnames(index2i,:));
      if s(1)=='_'
         s=['FIX',s]; % Matlab-struct names cannot start with _
      end
      indi=find(s=='(');
      for indii=indi'
         s(indii)='_';
      end;
      s(s==')'|s==' '|s=='['|s==']'|s==',')=[];
      if first
         data=struct(s,{val});
      else
         data=setfield(data,s,val);
      end
      i=i+length(index2);
      first=0;
   end
end

function data = unpackElement(dymstr,index)
% Scalar case.
% By Matthis Langelaar.
dataInfo=dymstr.dataInfo(index,:);
datamatrix=dataInfo(1); if datamatrix==0, datamatrix=2; end; % fix for time
datacol=abs(dataInfo(2));
datasign=sign(dataInfo(2));
% dataInfo(3) and (4) are ignored
eval(['data=datasign*dymstr.data_',int2str(datamatrix),'(:,datacol);']);

function [index,typ]=myStrMatch(name,strs)
% Written by Hans Olsson Dynasim.
% strmatch is not sufficient since
% a can match alpha, a.b, a[1]
%
% typ is:
%  0 scalar variable
%  1 array
%  2 model (struct)
index=strmatch(name,strs);
if length(name)==size(strs,2)
   typ=0;
   return;
end
indexOut=zeros(length(index),1);
j=0;
typ=0;
parlevel=0;
for i=1:length(index)
   c=strs(index(i),length(name)+1);
   if c=='('
      parlevel=parlevel+1;
   elseif c==')'
         parlevel=parlevel-1;
   elseif parlevel>0 | isletter(c) | c==',' | c=='_' | (c>='0' & c<='9')
      ;
   else
      j=j+1;
      indexOut(j)=index(i);
      if c=='[' 
         typ=1;
      elseif c=='.'
         typ=2;
      end
   end
end
index=indexOut(1:j);