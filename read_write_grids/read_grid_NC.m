function sData = read_grid_NC(path)


% varRem = {'lon', 'longitude', 'lat', 'latitude', 'time', 'time_bnds'};
% vatLat = {'lat', 'latitude'};
% vatLon = {'lon', 'longitude'};

% %Find possible data variables:
% if isempty(varargin) && isempty(varargin{1})
%     fInfo = ncinfo(path);
%     varLd = extractfield(fInfo, 'Variables');
%     if iscell(varLd) && isstruct(varLd{1})
%        varLd = extractfield(varLd{1},'Name'); 
%     end
% 
%     for ii = 1 : numel(varRem)
%         varLd(strcmpi(varLd, varRem{ii})) = [];
%     end
% else
%     vadLd = varargin{1};
% end

%Load variables:
fInfo = ncinfo(path);
varLd = extract_field(fInfo.Variables, 'Name');
% varLd = extractfield(fInfo.Variables, 'Name');

for ii = 1 : numel(varLd)
    sData.(varLd{ii}) = ncread(path, varLd{ii});
    
    attData = ncinfo(path, varLd{ii});
    if ~isempty(attData.Attributes)
        attData = squeeze(struct2cell(attData.Attributes))';
        sData.(['att' varLd{ii} ]) = attData;
        
        [rNoVal, ~] = ind2sub(size(attData),find(strcmpi(attData, 'missing_value') == 1));
        if ~isempty(rNoVal)
            valNoData = attData{rNoVal,2};
            sData.(varLd{ii})(sData.(varLd{ii}) == valNoData) = nan;
        end
        
        %Create date vector for time:
        if strcmpi(varLd{ii}, 'time')
            %Read time units from GCM file:
            [gcmRef, tUnits] = NC_time_units(attData);


            if regexpbl(gcmRef,'unknown')
                sData.date = nan(1,3);
            else
                %Find calendar info:
                strCal =  NC_cal(attData);

                if regexpbl(strCal, 'unknown')
                    strCal = 'gregorian';
                    warning('NC_avail_time:unknownCal','The calendar is unknown, but is being set to Gregorian.');
                end
                %Find time offset from the reference date:
                if regexpbl(tUnits, 'days since')
                    sData.date = days_2_date(sData.time, gcmRef, strCal);
                elseif regexpbl(tUnits, 'Ka BP') && regexpbl(strCal,'noleap')
                    sData.date = time_kaBP_2_standard(sData.time,gcmRef,strCal);
                elseif regexpbl(tUnits, 'hours since')
                    if numel(gcmRef) == 3
                        gcmRef = [gcmRef, 0];
                    end
                    sData.date = days_2_date(sData.time/24, gcmRef, strCal);
                else
                    disp(['The GCM' char(39) 's time units are' tUntInfo '.']);
                    error('NC_time:unitsUnknown','No case has been written to deal with the current time units.');
                end
            end
        end
    end
end

% for jj = 1 : varLon
%     if regexpbl(varLon{jj}, varLd)
%         
%     end
% end
% lonMODIS = ncread(path, 'lon');
%     lonMODIS = lonMODIS(:)';
% latMODIS = ncread(path, 'lat');
%     latMODIS = latMODIS(:);
% MODISTemp = ncread(path, varLd);

