from itertools import chain
import numpy as np
import xarray as xr

missing_gage = b'               '

@xr.register_dataset_accessor("routelink")
class Routelink:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        # Check that this is actually a routelink file? dimensions?

    def trace_links(
        self,
        index_or_gage,
        direction,
        max_depth=None
    ):

        if isinstance(index_or_gage, str):
            fw_gage = '{0: >{width}}'.format(index_or_gage, width=15).encode()
            the_index = np.where(self._obj['gages'] == fw_gage)[0].tolist()
        elif isinstance(index_or_gage, int):
            the_index = [index_or_gage]
        else:
            raise ValueError(
                'The index_or_gage argument is neither a str nor an int.')

        if direction is 'up':
            # Subtract 1 for zero-based/python indexing
            from_inds = self._obj.fromIndices.values - 1
            from_inds_start = self._obj.fromIndsStart.values - 1
            from_inds_end = self._obj.fromIndsEnd.values - 1
            def get_next_ups(index):
                return from_inds[
                    from_inds_start[index]:
                    from_inds_end[index]+1].tolist()
            get_next_trace = get_next_ups

        elif direction is 'down':
            # Subtract 1 for zero-based/python indexing
            to_ind = self._obj.toIndex.values - 1
            def get_next_down(index):
                next_down = [to_ind[index].tolist()]
                if next_down[0] is -1:
                    next_down = []
                return next_down
            get_next_trace = get_next_down

        else:
            raise ValueError("Direction must be one of 'up' or 'down'.")

        all_trace_inds = []
        trace_inds = the_index
        depth = 0

        while len(trace_inds) > 0:
            all_trace_inds += trace_inds
            if max_depth is not None and depth >= max_depth:
                break
            trace_inds = [get_next_trace(ii) for ii in trace_inds]
            trace_inds = list(chain.from_iterable(trace_inds))
            depth += 1

        _ = all_trace_inds.remove(the_index[0])
        return (the_index, all_trace_inds)

    def get_upstream_inds(self, *args, **kwargs):
        #index_or_gage, max_depth=None, get_gages=False, max_gages=None):
        return self.trace_links(direction='up', *args, **kwargs)
    #index_or_gage, direction='up', max_depth=max_depth,
    #        get_gages=False, max_gages=None)

    def get_downstream_inds(self, *args, **kwargs):
        #index_or_gage, max_depth=None, get_gages=False, max_gages=None):
        return self.trace_links(direction='down', *args, **kwargs)
    #index_or_gage, direction='down', max_depth=max_depth)
    
    def get_gages_from_inds(self, inds, keep_vars=['gages', 'feature_id']):
        rl_sub = self._obj.drop_vars(set(self._obj.variables) - set(keep_vars))
        rl_sub['index'] = xr.DataArray(
            np.arange(len(self._obj[keep_vars[0]])), dims='feature_id')
        rl_sub = rl_sub.isel({'feature_id': inds})
        rl_sub_2 = rl_sub.where(rl_sub.gages != missing_gage, drop=True)
        # This is annoying that where changes the types to float/object
        for var in list(rl_sub.variables):
            rl_sub_2[var] = rl_sub_2[var].astype(rl_sub[var].dtype)
        return rl_sub_2

