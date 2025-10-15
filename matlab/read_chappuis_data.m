function sigma_o3 = read_chappuis_data(wavelength_nm, T)
data_folder = '../data';
raw_data = readtable(fullfile(data_folder, 'SerdyuchenkoGorshelevVersionJuly2013-clip.csv'));

raw_T = [293, 283, 273, 263, 253, 243, 233, 223, 213, 203, 193];
raw_lam_nm = raw_data.vac;
sigma_o3 = interp2(raw_T, raw_lam_nm, raw_data{:,2:end}, T(:)', wavelength_nm(:), 'linear', 0);
sigma_o3 = sigma_o3 * 1e-4; % from cm^2 to m^2
end