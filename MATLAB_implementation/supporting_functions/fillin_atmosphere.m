% ATMOSPHERIC MODEL (pre-computed JB2008)
% valid from March 2008 - Feb 2224 for density file 'dens_jb2008_032020_022224.mat'
density_profile = 'JB2008'; % The options are 'static' or 'JB2008'
cfgMC.density_profile = density_profile;
if strcmpi(cfgMC.density_profile,'JB2008')
    cfgMC = initJB2008(cfgMC);
end

function cfgMCout = initJB2008(cfgMC)
    
    % Try multiple ways to find the density file
    fn = which('dens_jb2008_032020_022224.mat');
    if isempty(fn)
        % If which() doesn't find it, try relative paths
        possible_paths = {
            'dens_jb2008_032020_022224.mat',
            '../supporting_data/dens_jb2008_032020_022224.mat',
            '../../supporting_data/dens_jb2008_032020_022224.mat',
            '../../../supporting_data/dens_jb2008_032020_022224.mat'
        };
        
        for i = 1:length(possible_paths)
            if exist(possible_paths{i}, 'file')
                fn = possible_paths{i};
                break;
            end
        end
        
        if isempty(fn)
            error('Could not find dens_jb2008_032020_022224.mat. Please check the file exists in supporting_data directory.');
        end
    end
    
    load(fn);
    
    dens_times=zeros(length(dens_highvar.month),1);
    
    for k=1:length(dens_highvar.month)
        dens_times(k,1)= juliandate(datetime(dens_highvar.year(k),dens_highvar.month(k),0));
    end
    
    [dens_times2,alt2] = meshgrid(dens_times,dens_highvar.alt);

    cfgMCout = cfgMC;
    cfgMCout.param.dens_times = dens_times2;
    cfgMCout.param.alt = alt2;
    cfgMCout.param.dens_value = dens_highvar.dens;
end