function reduce_all_netcdfs_in_folder_to_checksum( foldername)

foldername = 'UPSY/component_test_mesh_creation/results';

henk = dir( foldername);

for i = 1: length( henk)
  filename = [foldername '/' henk(i).name];
  if endsWith( filename,'.nc')
    reduce_netcdf_to_checksum( filename);
    delete( filename)
  end
end

end