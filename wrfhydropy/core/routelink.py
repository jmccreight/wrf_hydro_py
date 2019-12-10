from itertools import chain
import numpy as np
import xarray as xr

@xr.register_dataset_accessor("routelink")
class Routelink:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._center = None

    def get_upstream_links(self, index_or_gage):
        if isinstance(index_or_gage, str):
            fw_gage = '{0: >{width}}'.format(index_or_gage, width=15).encode()
            the_index = np.where(self._obj['gages'] == fw_gage)[0].tolist()
        elif isinstance(index_or_gage, int):
            the_index = [index_or_gage]
        else:
            raise ValueError(
                'The index_or_gage argument is neither a str nor an int.')

        # Subtract 1 for zero-based/python indexing
        # to_inds = self._obj.toIndex.values - 1
        from_inds = self._obj.fromIndices.values - 1
        from_inds_start = self._obj.fromIndsStart.values - 1
        from_inds_end = self._obj.fromIndsEnd.values - 1

        def get_next_ups(index):
            return from_inds[
                from_inds_start[index]:
                from_inds_end[index]+1].tolist()

        all_upstream_inds = []
        up_inds = the_index
        while len(up_inds) > 0:
            all_upstream_inds += up_inds
            up_inds = [get_next_ups(ii) for ii in up_inds]
            up_inds = list(chain.from_iterable(up_inds))

        _ = all_upstream_inds.remove(the_index[0])
        return (the_index, all_upstream_inds)
