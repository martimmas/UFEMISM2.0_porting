function compare_netcdf( filename_ref, filename_mod)
% Compare two NetCDF files to see if they are identical,
% and if not, point out where the differences are

filename_ref_short = shorten_filename( filename_ref);
filename_mod_short = shorten_filename( filename_mod);
if ~strcmpi( filename_ref_short, filename_mod_short)
  error('files should have the same name')
end

disp(['Comparing ' filename_ref_short])

info_matches = compare_ncinfo( filename_ref, filename_mod);
if ~info_matches
  error('NetCDF info does not match')
else
  data_matches = compare_data( filename_ref, filename_mod);
end
if ~data_matches
  error('data in file does not match')
end

  function filename_short = shorten_filename( filename)
    i = strfind( filename,'/');
    if isempty(i)
      filename_short = filename;
    else
      filename_short = filename( i(end)+1:end);
    end
  end

  function info_matches = compare_ncinfo( filename_ref, filename_mod)
    f_ref = ncinfo( filename_ref);
    f_mod = ncinfo( filename_mod);
    nerrs = compare_ncinfo_structs( 0, shorten_filename( filename_ref), [], f_ref, f_mod);
    info_matches = nerrs == 0;
  end

  function nerrs = compare_ncinfo_structs( nerrs, basename, fieldname, s1, s2)

    % disp( basename)

    if ~strcmpi( class( s1), class( s1))
      nerrs = nerrs + 1;
      disp(['Mismatching classes in ' basename])
      return
    end

    if isnumeric( s1) || islogical( s1)
      eq = all( s1 == s2);
      if ~eq
        nerrs = nerrs + 1;
        disp(['Mismatch in ' basename])
      end
      return
    end
    
    if ischar( s1)
      eq = strcmpi( s1, s2);
      % Exceptions
      if strcmpi( fieldname,'Filename')
        eq = true;
      end
      if ~eq
        nerrs = nerrs + 1;
        disp(['Mismatch in ' basename])
      end
      return
    end

    if ~isstruct( s1)
      error(['invalid class "' class( s1) '"'])
    end

    % Now we know s1,s2 are structs

    fieldnames1 = fields( s1);
    fieldnames2 = fields( s2);

    if length( fieldnames1) ~= length( fieldnames2)
      nerrs = nerrs + 1;
      disp(['Mismatching number of fields in ' basename])
      return
    end

    for fi = 1: length( fieldnames1)
      fieldname1 = fieldnames1{ fi};
      fieldname2 = fieldnames2{ fi};
      if ~strcmpi( fieldname1, fieldname2)
        nerrs = nerrs + 1;
        disp(['Mismatching field names in ' basename])
        return
      end
    end

    for fi = 1: length( fieldnames1)

      fieldname = fieldnames1{ fi};

      field1 = s1.(fieldname);
      field2 = s2.(fieldname);

      if isstruct( field1) && length( field1) > 1
        for fii = 1: length( field1)
          nerrs = compare_ncinfo_structs( nerrs, [basename '/' fieldname], fieldname, ...
            field1( fii), field2( fii));
        end
      else
        nerrs = compare_ncinfo_structs( nerrs, [basename '/' fieldname], fieldname, ...
          field1, field2);
      end
  
    end

  end

  function data_matches = compare_data( filename_ref, filename_mod)

    data_matches = true;

    f = ncinfo( filename_ref);

    % Regular variables;
    for vi = 1: length( f.Variables)
      var_name = f.Variables( vi).Name;
      data_matches = data_matches && compare_variable( ...
        filename_ref, filename_mod, var_name);
    end

    % Grouped variables
    for gi = 1: length( f.Groups)
      group_name = f.Groups(gi).Name;
      for vi = 1: length( f.Groups(gi).Variables)
        var_name = [group_name '/' f.Groups(gi).Variables( vi).Name];
        data_matches = data_matches && compare_variable( ...
          filename_ref, filename_mod, var_name);
      end
    end
  end

  function data_matches = compare_variable( filename_ref, filename_mod, var_name)

    d_ref = ncread( filename_ref, var_name);
    d_mod = ncread( filename_mod, var_name);

    dd = all( d_ref == d_mod);
    data_matches = all( dd(:));

    if ~data_matches
      disp(['Mismatch in ' var_name])
    end

  end

end