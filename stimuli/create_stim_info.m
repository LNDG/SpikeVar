function info = create_stim_info(image_dir, save_dir)
%create_stim_info Creates table with information regarding the stimuli.
%   Takes as input the directory containing the stimuli

%% set up version, categroy names etc
versions = ["newolddelay", "newolddelay2", "newolddelay3"];
cat_dict = [{"houses", "landscapes", "mobility", "phones", "smallAnimal"},...
    {"fruit", "kids", "military", "space", "zzanimal"},...
    {"1cars", "2food","3people", "4spatial", "5animals"}];
cat_nums = 1:15;

%% loop across versions
stim_codes = []; pic_paths = [];
stim_versions = []; stim_categories = [];
for ver = 1:length(versions)
    verdir = [image_dir + versions(ver) + '/'];
    % loop across categories
    for cat = 1:5
        catdir = [verdir + cat_dict(5*(ver-1)+cat) + '/'];
        % get list of images in the folder of interest
        piclist = dir(catdir + '*.jpg');
        npics = length(piclist);
        % loop across pictures in folder, load and make grayscale
        for ipic = 1:npics
            % store version and category for each stimulus 
            stim_versions = [stim_versions; ver];
            stim_categories = [stim_categories; cat];
            % get the unique stimulus number that is used later
            stim_code = 1000 + cat_nums(5*(ver-1)+cat)*100+ipic;
            stim_codes = [stim_codes; stim_code];
            % load image
            pic_path = [catdir + piclist(ipic).name];
            pic_paths = [pic_paths; pic_path];
        end
    end
end
%% prepare info data and save
info = table(NaN(590,1), NaN(590,1), NaN(590,1), ... 
    'VariableNames', {'stimulus_code', 'version', 'category'});
info.stimulus_code = stim_codes;
info.version = stim_versions;
info.category = stim_categories;
info.picture_path = pic_paths;
save([save_dir + 'info.mat'], 'info')

end

