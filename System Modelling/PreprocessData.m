function ProcessedData = PreprocessData(data,TestLen,varargin)
    FieldNames = data.getElementNames;
    FieldNames = FieldNames(~cellfun('isempty',FieldNames));
    for ii = 1:length(FieldNames)
        foo = data.getElement(FieldNames{ii});
        foo = squeeze(foo.Values.Data);
        ProcessedData.(FieldNames{ii}) = foo;
    end
    t = data.getElement(FieldNames{1}); 
    t = t.Values.time;
    ProcessedData.('time') = t;
    
    StructNames = fieldnames(ProcessedData);
    for ii = 1:numel(StructNames)
        lens(ii) = length(ProcessedData.(StructNames{ii}));
    end
    for ii = 1:numel(StructNames)
        if length(ProcessedData.(StructNames{ii})) < max(lens)
            for jj = 1:size(ProcessedData.(StructNames{ii}),2)
                try
                ProcessedData.(StructNames{ii})(1:max(lens),jj) = interp1(1:max(lens)/length(ProcessedData.(StructNames{ii})):max(lens),ProcessedData.(StructNames{ii})(:,jj),1:1:max(lens),'previous','extrap')';
                catch
                sprintf('Issue with %s',StructNames{ii})
                end
            end
        end
    end
    if nargin == 2
        for ii = 1:numel(StructNames)
            for jj = 1:size(ProcessedData.(StructNames{ii}),2)
            ProcessedData.(StructNames{ii})(1:TestLen,jj) = ProcessedData.(StructNames{ii})(1:TestLen,jj);
            end
        end
    end
end