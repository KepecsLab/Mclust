function [T,WV] = HDF5LoadingEngine(fn, R, record_units)

persistent T0 WV0 fn0

if isempty(T0) || ~strcmp(fn, fn0)
	fn0 = fn;
	T0 = hdf5read(fn, '/t');
	units = hdf5read(fn, '/unit');
	wv = -hdf5read(fn, '/wv');
	WV0 = permute(wv, [3 2 1]);
end

switch(record_units)
	case 1 % records_to_get is timestamp list
		keep = false(size(T0))
		for iT = 1:length(R)
			f = find(R(iT)==T0,1,'first');
			keep(f) = true;
		end
		T = T0(keep);
		WV = WV0(keep,:,:);
	case 2 % record number list
		T = T0(R);
		WV = WV0(R,:,:);
	case 3 % records to get is range of timestamps
		tstart = find(T0 >= R(1), 1, 'first');
		tend = find(T0 <= R(end), 1, 'last');
		if ~isempty(tstart) & ~isempty(tend)
			T = T0(tstart:tend);
			WV = WV0(tstart:tend,:,:);
		else
			T = [];
			WV = [];
		end
	case 4 % records to get is range of records
		tstart = R(1); tend = R(end);
		T = T0(tstart:tend);
		WV = WV0(tstart:tend,:,:);
	case 5 % length
		T = length(T0);
		WV = [];
end