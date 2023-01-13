function age = get_ages()
    demo_dir = '/Users/kamp/PhD/spikevar/data/demographics/';
    age = load([demo_dir 'AgeSex.mat']).corval;
    age = sortrows(age,1);
end