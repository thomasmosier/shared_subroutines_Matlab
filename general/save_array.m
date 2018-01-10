function save_array(path, array, strArray)

eval([strArray '= array;']);
save(path, strArray);