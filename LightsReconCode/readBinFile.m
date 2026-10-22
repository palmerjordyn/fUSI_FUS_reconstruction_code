% function A = readBinFile(fname, datasize, datatype)
% 
%     if ~exist('datatype', 'var') || isempty(datatype)
%         datatype = 'double';
%     end
% 
%     fileID = fopen(fname);
%     if strcmp(datatype,'complex')
%         try
%             %A = freadcomplex(fileID, prod(datasize), 'double');
%             raw = fread(fileID, 2*n, precision);
%             A = complex(raw(1:2:end), raw(2:2:end));
%         catch ME
%             fclose(fileID);
%             rethrow(ME)
%         end
%     else
%         try
%             A = fread(fileID, prod(datasize), datatype);
%         catch ME
%             fclose(fileID);
%             rethrow(ME)
%         end
%     end
%     fclose(fileID);
%     A = reshape(A, datasize);
% end


function A = readBinFile(fname, datasize, datatype)

    if nargin < 3 || isempty(datatype)
        datatype = 'double';
    end

    fileID = fopen(fname, 'rb');

    if fileID == -1
        error('Could not open file: %s', fname);
    end

    cleanup = onCleanup(@() fclose(fileID));

    n = prod(datasize);

    if strcmp(datatype, 'complex')
        raw = fread(fileID, 2*n, 'double=>double');

        if numel(raw) ~= 2*n
            error('File does not contain the expected number of complex samples.');
        end

        A = complex(raw(1:2:end), raw(2:2:end));

    else
        A = fread(fileID, n, [datatype '=>' datatype]);

        if numel(A) ~= n
            error('File does not contain the expected number of samples.');
        end
    end

    A = reshape(A, datasize);

end