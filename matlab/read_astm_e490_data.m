function E = read_astm_e490_data(wavelength_nm)
data_folder = '../data';
raw_data = readtable(fullfile(data_folder, 'e490_00a_amo.csv'), 'VariableNamingRule', 'preserve');
E = interp1(raw_data{:,1}, raw_data{:,2}, wavelength_nm*1e-3, 'linear', 0);  % W/m^2/um
E = E(:)*1e-3;  % convert to W/m^2/nm
end