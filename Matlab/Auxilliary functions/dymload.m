function [dymstr,err] = dymload(matfile)
% dymload  load dymola simulation result into Matlab workspace
% 
%    d = dymload;        loads data from 'dsres.mat'.
%    d = dymload(file);  loads data from '<file>.mat'.
%
%    d is a Matlab STRUCT-variable containing the following fields:
%    .fname        the filename
%    .pname        the pathname
%    .nnames       the number of variable names
%    .ndatamat     the number of data matrices
%    .name         the names of the variables (each row a name)
%    .description  the descriptions of the variables (each row a description)
%    .dataInfo     the dataInfo-array (see Dymola output format specification)
%    .data_#       the data-matrices. # ranges from 1 to .ndatamat
%
%    Only the Dymola simulation result data format version 1.1 is supported,
%    that is compatible with Dymola 3.1 and higher.
%
% See also: dymbrowse, dymget

% # this function is a slightly altered version of the tload-function
% Copyright (c) 2000 DLR,
% Copyright (c) 2001 Dynasim AB, Sweden
%    All rights reserved.

% determine file name
err = '';
dymstr=[];
  if nargin < 1
     file = 'dsres.mat';
  else
     [f_p, f_n, f_e]=fileparts(matfile);
     if isempty(f_e)
       file=[matfile,'.mat'];
     elseif strcmp(lower(f_e),'.mat') 
       file=matfile;
     else
       err = ['filename (= "',matfile,'") has not extension ".mat"'];
        if nargout<2
           error(err);
        else
           return;
        end
     end
  end

% check file existence
   if ~size(dir(file),1)
      error( sprintf( ['"', file, '" does not exist in \n   %s'], pwd) )
   end

% load data
  load(file);

% read Aclass variable
  matlabVersion = version;
  if exist('Aclass') ~= 1
     if matlabVersion(1,1)=='4'
        if exist('class') ~= 1
           err = ['no traj. on file "' file '" ("Aclass" is missing).'];
           if nargout<2
              error(err);
           else
              return;
           end
        else
           Aclass = class;
        end
     else
        err =['no trajectory on file "' file '" ("Aclass" is missing).'];
        if nargout<2
          error(err);
        else
          return;
        end
     end
  end

% check whether file has correct class name
  classReq = 'Atrajectory';
  ncol1 = size(classReq,2);
  [nrow2,ncol2] = size(Aclass);
  if ncol1 < ncol2 
     classReq = [ classReq, blanks(ncol2-ncol1) ];
  elseif ncol1 > ncol2
     Aclass = [ Aclass, blanks(ncol1-ncol2) ];
     ncol2  = size(Aclass, 2);
  end
  if nrow2 < 2 then
     err = [ 'file "' file '" is not of class ' classReq ];
     if nargout<2
       error(err);
     else
       return;
     end

  elseif Aclass(1,:) ~= classReq(1,:)
     err = [ 'file "' file '" is not of class ' classReq ];
     if nargout<2
       error(err);
     else
       return;
     end
  end

% Check version number
  if ['1.0'] == Aclass(2,1:3)
     vers = 0;
  elseif ['1.1'] == Aclass(2,1:3)
     vers = 1;
  else
     err = [ 'file "' file '" has wrong version number ' Aclass(2,:) ];
     if nargout<2
       error(err);
     else
       return;
     end
  end

% Determine whether matrices have to be transposed
  if nrow2 < 4
     trans = 0;
  elseif Aclass(4,1:8) == ['binTrans']
     trans = 1;
  else
     trans = 0;
  end
  
  out.fname = file;
  out.pname = pwd;
  
% Action according to version number
  if vers == 0 
     % Refuse (not implemented because no data format information available)
     % error( ['Wrong version number (version 1.1 required)'] )
     if exist('names') ~= 1
        err = ['no traj. on file "' file '" (matrix "names" missing).'];
        if nargout<2
          error(err);
        else
          return;
        end
     end
     if trans == 0
        out.name = names;
     else
        out.name = names';
     end;
     out.nnames=size(out.name,1);
     out.ndatamat=1;
     if trans == 0
        out.data_1 = data;
     else
        out.data_1 = data';
     end;
     out.description = blanks(out.nnames)';
     out.dataInfo = [ones(out.nnames,1),(1:out.nnames)',zeros(out.nnames,1),-1*ones(out.nnames,1)];
  else
     % Check existance of name and dataInfo matrix
       if exist('name') ~= 1
          err =['no traj. on file "' file '" (matrix "name" xmissing).'];
          if nargout<2
            error(err);
          else
            return;
          end
       elseif exist('dataInfo') ~= 1
          err = ['no traj. on file "' file '" ("dataInfo" missing).'];
          if nargout<2
            error(err);
          else
            return;
          end
       end

     % Copy name
       if trans == 0
          out.name = name;       
          out.dataInfo = dataInfo;
       else
          out.name = name';
          out.dataInfo = dataInfo';
       end
       if exist('description') ~=1
          out.description=blanks(out.nnames)';
       elseif trans == 0
          out.description = description;   
       else
          out.description = description';
       end
       out.nnames=size(out.name,1);
       out.ndatamat=max(out.dataInfo(:,1));
       
     % Store matrices data_i
       for i=1:out.ndatamat
          % Determine matrix name
            if trans == 0
               eval( ['out.data_' int2str(i) '= data_' int2str(i) ';'] );
            else
               eval( ['out.data_' int2str(i) '= data_' int2str(i) ''';'] );
            end
       end
  end

% print info message
  if nargout>=2
   
  elseif file(1,1) == '/' | file(1,1) == '\' | file(1,2) == ':'
     disp( [ '> ', file, ' loaded.' ] )
  else
     machine = computer;
     if machine(1,1:2) == 'PC'
        disp( [ '> ', lower(pwd), '\', lower(file), ' loaded.'] )
     else
        disp( [ '> ', pwd, '/', file, ' loaded.' ] )
     end
  end
  
  dymstr=out;