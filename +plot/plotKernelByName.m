function plotKernelByName(name, fitdirbase, outdir)
% dts = {'20121104-LIP-1', '20130611-LIP-1', '20140307-LIP-4', ...
%     '20150331-LIP-5', '20150331-LIP-7', '20150401-LIP-10', ...
%     '20140218-MT-1', '20140304-MT-2', '20150324a-MT-10', ...
%     '20150324a-MT-13'};

    if nargin < 3
        outdir = 'tmp';
    end
    if nargin < 2
%         fitdirbase = '20150615';
        fitdirbase = '20150716';
%         fitdirbase = 'sametargs';
    end
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    tags = strsplit(name, '-');
    dt = tags{1};
    lbl = strjoin(tags(2:end), '_');
    isNancy = str2num(dt(4)) > 4;

    if isNancy
        mnkNm = 'nancy';
    else
        mnkNm = 'pat';
    end
    fitdir = ['data/' fitdirbase '-' mnkNm '/fits'];

    d = io.loadDataByDate(dt, isNancy);
    fs = io.loadFitsByDate(dt, fitdir);
    f = fs.(lbl).ASD{end}; 
    plot.plotAndSaveKernel(f, d, outdir, true, true, false);
end
