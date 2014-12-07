% Code to determine structure of an SPM.mat file.

% An SPM object must be a 1x1 struct.  SPM accomplishes that by placing
% each item that corresponds to a field within a single item cell array.

% e.g.: struct('a',{struct('a',6,'b',7)},'b',{{1,'horse','banana'}},'c',{9})

function spm_iterator(SPM, parent)
    
    if ~ischar(parent)
        error('Argument ''parent'' must be a character array');
    end
    
    if ~isstruct(SPM)
        error('Argument ''SPM'' is not a structural array.');
    end
    fields = fieldnames(SPM);

    for i=1:length(fields)
      disp(strcat(parent, '.', fields{i}, ' : ', class(SPM.(fields{i}))));
      if isstruct(SPM.(fields{i})) && size(SPM.(fields{i}),2)>1
          disp('Structural array field names: ');
          for field=fieldnames(SPM.(fields{i}))
              disp(field);
          end
      end
      if iscell(SPM.(fields{i}))
          disp('Types of each cell element: ');
          tmp=SPM.(fields{i});
          for elm=1:length(tmp)
              disp(class(tmp{elm}));
          end
      end
      disp(' ');
      if isstruct(SPM.(fields{i})) && size(SPM.(fields{i}),1)==1 && size(SPM.(fields{i}),2)==1
          spm_iterator(SPM.(fields{i}), strcat(parent, '.', fields{i}));
      end
    end
end