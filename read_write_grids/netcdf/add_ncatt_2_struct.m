function sGeo = add_ncatt_2_struct(sGeo, pathLd, varCurr)
    [attCurr, valMiss] = load_ncatt(pathLd, varCurr);
    if ~isempty(attCurr)
        sGeo.(['att' varCurr]) = attCurr;
    end
    
    if ~isnan(valMiss)
       sGeo.(varCurr)( sGeo.(varCurr) == valMiss ) = nan;  
    end
end