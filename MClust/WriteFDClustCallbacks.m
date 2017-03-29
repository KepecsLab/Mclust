function WriteFDClustCallbacks

% CBP_Callbacks
%
% 

    cboHandle = gcbo;
    figHandle = gcf;
    callbackTag = get(cboHandle, 'Tag');
  
global MClust_FeatureNames MClust_FeatureSources MClust_FDfn MClust_FeatureTimestamps
global MClust_Clusters MClust_Colors

% get parameters from figure
newFeatures = get(findobj(figHandle,'Tag','FeaturesUseListbox','TooltipString', ...
    'These features will be written out.'),'String');

iClust = get(findobj('Tag','WriteFeatureData'), 'UserData');
close(figHandle);

% get features, construct FeatureData
[spikeIndex MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});

if isempty(newFeatures)
    errordlg('No features chosen', 'WriteFDData Error', 'modal');
    return;
else
	nFeatures = length(newFeatures);
	FeatureData = nan(length(spikeIndex), nFeatures);
	for iF = 1:nFeatures
		FeatureID = strmatch(newFeatures{iF}, MClust_FeatureNames);
		if isempty(FeatureID)
			errordlg('No matching feature found', 'WriteFDData Error', 'modal');
			return;
		elseif length(FeatureID) > 1
			errordlg('Too many matching features found', 'WriteFDData Error', 'modal');
			return;
		else
			temp = load(MClust_FeatureSources{FeatureID,1}, '-mat', 'FeatureData');
			FeatureData(:,iF) = temp.FeatureData(spikeIndex,MClust_FeatureSources{FeatureID,2});
		end
	end
end

[fpath fname fext] = fileparts(MClust_FDfn);
FETfn = fullfile(fpath,[fname '-Cluster' num2str(iClust) '.mat']);
TS = MClust_FeatureTimestamps(spikeIndex);
save(FETfn, 'FeatureData', 'TS');
