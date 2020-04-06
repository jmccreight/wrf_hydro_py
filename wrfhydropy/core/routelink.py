from itertools import chain
import numpy as np
import xarray as xr

#                123456789012345
missing_gage = b'               '

@xr.register_dataset_accessor("routelink")
class Routelink:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        # Check that this is actually a routelink file? dimensions?

    def trace_link(
        self,
        index_or_gage,
        direction: str,
        max_depth: int = None,
        id_in: bool = False
    ):

        if isinstance(index_or_gage, str):
            fw_gage = '{0: >{width}}'.format(index_or_gage, width=15).encode()
            the_index = np.where(self._obj['gages'] == fw_gage)[0].tolist()
        elif isinstance(index_or_gage, (int, np.integer)):
            if id_in:
                the_index = np.where(self._obj['link'] == index_or_gage)[0].tolist()
            else:
                the_index = [index_or_gage]
            if the_index[0] >= len(self._obj.link):
                raise ValueError(
                    "You appear to have passed a id and not an index. This"
                    "requires passing the argument 'id=True'.")
        else:
            raise ValueError(
                'The index_or_gage argument is neither a str nor an int.')

        if direction == 'up':
            # Subtract 1 for zero-based/python indexing
            from_inds = self._obj.fromIndices.values - 1
            from_inds_start = self._obj.fromIndsStart.values - 1
            from_inds_end = self._obj.fromIndsEnd.values - 1
            def get_next_ups(index):
                return from_inds[
                    from_inds_start[index]:
                    from_inds_end[index]+1].tolist()
            get_next_trace = get_next_ups

        elif direction == 'down':
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

    def trace_links(
        self,
        index_or_gage,
        direction: str,
        max_depth: int = None,
        id_in: bool = False
    ):
        if isinstance(index_or_gage, (str, int, np.integer)):
            return self.trace_link(
                index_or_gage, direction, max_depth=max_depth, id_in=id_in)
        elif isinstance(index_or_gage, xr.DataArray):
            index_or_gage = index_or_gage.values.tolist()
        elif isinstance(index_or_gage, xr.DataArray):
            index_or_gage = index_or_gage.values.tolist()
        return [
            self.trace_link(ii, direction, max_depth=max_depth, id_in=id_in)
            for ii in index_or_gage]

    def get_upstream_inds(self, *args, **kwargs):
        return self.trace_links(direction='up', *args, **kwargs)

    def get_upstream_ids(self, *args, **kwargs):
        return self.inds_to_ids(
            self.trace_links(direction='up', *args, **kwargs)[1])

    def get_downstream_inds(self, *args, **kwargs):
        return self.trace_links(direction='down', *args, **kwargs)

    def get_downstream_ids(self, *args, **kwargs):
        inds_result = self.trace_links(direction='down', *args, **kwargs)
        return [(self.inds_to_ids(ii[0]), self.inds_to_ids(ii[1]))
                for ii in inds_result]

    def id_to_ind(self, id_in: int):
        if not isinstance(id_in, (int, np.integer)):
            raise ValueError('Input argument must be integer.')
        return np.where(self._obj['link'] == id_in)[0][0]

    def ids_to_inds(self, id_in: list):
        if isinstance(id_in, (int, np.integer)):
            id_in = [id_in]
        elif not isinstance(id_in, list):
            raise ValueError('Input argument must be list or integer.')
        return [self.id_to_ind(ii) for ii in id_in]

    def ind_to_id(self, ind_in: int):
        if not isinstance(ind_in, (int, np.integer)):
            raise ValueError('Input argument must be integer.')
        return self._obj['link'].isel(feature_id=ind_in).values.tolist()

    def inds_to_ids(self, ind_in: list):
        if isinstance(ind_in, (int, np.integer)):
            ind_in = list(ind_in)
        elif not isinstance(ind_in, list):
            raise ValueError('Input argument must be list or integer.')
        return [self.ind_to_id(ii) for ii in ind_in]

    def get_outlet_inds(self, index_or_gage, id_in=False):
        # this could use a generator since only the last value is kept
        if isinstance(index_or_gage, (int, np.integer, str)):
            index_or_gage = [index_or_gage]
        elif not isinstance(index_or_gage, list):
            raise ValueError('Input argument must be list or integer.')       
        inds = self.trace_links(index_or_gage, id_in=id_in, direction='down')
        outlets = []
        for ii in range(len(inds)):
            # if it is its own outlet!
            if len(inds[ii][1]) == 0:
                outlets += [(inds[ii][0], inds[ii][0])]
            else:
                outlets += [(inds[ii][0], [inds[ii][1][-1]])]
        return outlets

    def get_outlet_ids(self, index_or_gage, id_in=False):
        inds = self.get_outlet_inds(index_or_gage, id_in=id_in)
        outlets = []
        for ii in range(len(inds)):
            outlets += [
                (self.inds_to_ids(inds[ii][0]),
                 self.inds_to_ids(inds[ii][1]))]
        return outlets

    def id_to_gage(self, id_in: int):
        return self.ind_to_gage(self.id_to_ind(id_in))

    def ids_to_gages(self, ids_in: list, drop_missing: bool = True):
        ret_val = self.inds_to_gages(
            self.ids_to_inds(ids_in), drop_missing=drop_missing)
        if drop_missing:
            ret_ids = self.inds_to_ids(ret_val[0])
            return (ret_ids, ret_val[1])
        else:
            return ret_val

    def ind_to_gage(self, ind_in: int):
        if not isinstance(ind_in, (int, np.integer)):
            raise ValueError('Input argument must be integer.')
        return self._obj['gages'].isel(feature_id=ind_in).values.tolist()

    def inds_to_gages(self, ind_in: list = [], drop_missing: bool = True):
        if ind_in == []:
            ind_in = rl._obj.feature_id.values.tolist()
        if isinstance(ind_in, (int, np.integer)):
            ind_in = list(ind_in)
        elif not isinstance(ind_in, list):
            raise ValueError('Input argument must be list or integer.')
        gage_list = [self.ind_to_gage(ii) for ii in ind_in]
        if drop_missing:
            gage_list = [gg for gg in gage_list if gg != missing_gage]
            gage_inds = self._obj.feature_id.where(
                self._obj.gages.isin(gage_list), drop=True).values.tolist()
            gage_inds = [int(gg) for gg in gage_inds]
            return (gage_inds, gage_list)
        else:
            return gage_list

    def gage_inds_by_outlet_ind(gage_inds: list = []):
        if gage_inds == []:
            gage_inds = self.inds_to_gages()
        outlet_gages = {}
        for oo in outlets:
            outlet_gages[oo] = []
            for go in gage_outlets:
                if go[1][0] == oo:
                    outlet_gages[oo] += [go[0][0]]
        
    
    def get_nested_gages():

        







# # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# # Copyright UCAR (c) 2018
# # University Corporation for Atmospheric Research(UCAR)
# # National Center for Atmospheric Research(NCAR)
# # Research Applications Laboratory(RAL)
# # P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# # 17/07/2018
# #
# # Name:        module1
# # Purpose:
# # Author:      $ Kevin Sampson
# # Created:     17/07/2018
# # Licence:     <your licence>
# # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

# '''
# 7/17/2018:
#     This script will retroactively add the 'From' vector and mask vector to an
#     existing RouteLink file in order to allow the upstream/downstream network
#     tracing to be performed more effiently within NMW/WRF-Hydro code.

# 7/18/2018:
#     The script was modified to produce a set of vectors that describe the from
#     flowlines for each link in an index-based fashion. This method uses a vector
#     to describe the start and stop index for finding all upstream flowlines for
#     any individual flowline.

# '''

# # Import Modules
# import sys
# import os
# import time
# import numpy
# import netCDF4
# from collections import OrderedDict

# # Globals

# # File inputs
# #inNC = r'C:\Data\Projects\Gochis\NWM_v2_0\In_Use\RouteLink_NWMv2.0_FromVars.nc'      # Input file to add ascendingIndex to
# dim = 'feature_id'                                                                 # Dimension name for input variable
# IDvar = 'link'                                                                  # Variable name containing the ID to be indexed over
# TOvar = 'to'                                                                    # Variable indicating downstream topology

# # Outputs
# fromVar = 'fromIndices'                                                         # Variable name to store the from segment list
# toIndexVar = 'toIndex'

# # Options
# addDim = True                                                                   # Switch for adding a new dimension
# removeZeros = True                                                              # Switch to remove 0 as a legitimate ID value. If True, addDim must be True also
# maskStyle = False                                                               # Switch for which method of defining upstream segments. This is the mask-vector method
# indexStyle = True                                                               # Switch for which method of defining upstream segments. This is the indexed array method

# # Functions
# def checkIndices(rootgrp, Tos):

#     tic1 = time.time()

#     # Obtain all input arrays
#     idsVar = rootgrp.variables[IDvar][:]
#     toVar = rootgrp.variables[TOvar][:]
#     FromVar = rootgrp.variables[fromVar][:]
#     StartVar = rootgrp.variables['fromIndsStart'][:]
#     EndVar = rootgrp.variables['fromIndsEnd'][:]
#     ToVar = rootgrp.variables[toIndexVar][:]

#     # Check the number of 'to' segment ComIDs that match between the toIndex variable and the existing 'to' variable
#     match1 = (idsVar[ToVar-1] == toVar).sum()
#     if match1 == toVar[toVar!=0].shape:
#         print('Sorted and indexed "link" array matches "to" array: %s/%s' %(match1, toVar[toVar!=0].shape[0]))

#     '''Test 1: Make sure that the toIndices relate back to the appropriate link IDs.'''
#     if (numpy.where(numpy.in1d(toVar, idsVar))[0] == numpy.where(numpy.in1d(toVar, idsVar[ToVar[ToVar!=0]-1]))[0]).sum() == FromVar.shape[0]:
#         print('The resorted IDs according to toIndices matches the existing "link":"to" variables.')

#     '''Test 2: Ensure that as you traverse through the Ids in "link", that the
#     indices make sense for the upstream and downstream flowlines.'''
#     samples = idsVar.shape[0]                                                   #samples = 5000
#     for Ind in range(samples):
#         linkID = idsVar[Ind]
#         upStart = StartVar[Ind]-1                                               # Convert 1-based index to 0-based
#         upEnd = EndVar[Ind]-1                                                   # Convert 1-based index to 0-based
#         downSeg = ToVar[Ind]-1 if ToVar[Ind]>0 else ToVar[Ind]                  # Convert 1-based index to 0-based except for values of 0
#         ups = idsVar[FromVar[upStart:upEnd+1]-1]
#         linkTO = idsVar[downSeg] if downSeg!=0 else 0
#         for item in ups:
#             if Tos.get(item) != linkID:
#                 print('[%s] Up fail: %s, up: %s' %(Ind, linkID, item))
#         if Tos.get(linkID) != linkTO:
#             print('[%s] Down fail: %s, %s' %(Ind, linkID, linkTO))
#     return

# if __name__ == '__main__':

#     tic = time.time()

#     # Read input file
#     inNC = sys.argv[1]
#     print('    Starting to read input RouteLink file: %s' %inNC)
#     rootgrp = netCDF4.Dataset(inNC, 'r+')
#     outDtype = rootgrp.variables[IDvar].dtype
#     origDim = rootgrp.variables[IDvar].dimensions
#     ids = rootgrp.variables[IDvar][:].tolist()
#     tos = rootgrp.variables[TOvar][:].tolist()
#     listlen = len(tos)
#     idsArr = numpy.array(ids)                                               # Make an array out of the link IDs

#     # Create from dictionary from To dictionary (for compatibility with existing RouteLink generation scripts)
#     Tos = OrderedDict([(fromID,toID) for fromID, toID in zip(ids[:], tos[:])])  # Build an ordered dictionary of downstream topology - slow
#     del ids, tos

#     # Invert the dictionary (fast). Make a dictionary of the upstream segments for each segment
#     inv_map = OrderedDict()
#     for k, v in Tos.items():
#         inv_map[v] = inv_map.get(v, [])
#         inv_map[v].append(k)

#     if removeZeros:
#         inv_map.pop(0)                                                          # Remove values of 0 (network endpoints)

#     # Flatten the lists
#     FromList = [item for sublist in inv_map.values() for item in sublist]       # The list containing the source flowline (unique)
#     FromwhichList = [key for key,val in inv_map.items() for item in val]    # The list containing the destination flowline (non-unique)
#     del inv_map

#     # Use the mask array method (similar to spatial weight files)
#     if maskStyle:
#         '''
#         You could either use the existing dimension and allow values of 0 in the fromwhich
#         variable, or remove these network endpoint values and shorten the length of these
#         new variables, but at the cost of adding a dimension.
#         '''

#         # Create new dimension and variables
#         fromMask = 'fromwhich'
#         newDim = 'fromDim'
#         if addDim:
#             dimName = newDim
#             rootgrp.createDimension(newDim, len(FromList))                      # Create Dimension
#         else:
#             dimName = origDim

#         From = rootgrp.createVariable(fromVar, outDtype, dimName)               # Variable (32-bit signed integer)
#         Fromwhich = rootgrp.createVariable(fromMask, outDtype, dimName)         # Variable (32-bit signed integer)
#         From.long_name = 'Vector of upstream segments for segments listed in fromwhich variable'
#         Fromwhich.long_name = 'Mask of flowline IDs for all upstream reaches listed in from variable'

#         # Populate newly created variables
#         From[:] = numpy.array(FromList)
#         Fromwhich[:] = numpy.array(FromwhichList)

#     # Use the index start:stop style, like MizuRoute
#     elif indexStyle:
#         '''
#         This block will use the from and to information to build an array of all
#         from and to flowlines as two separate arrays (FromArr, FromwhichArr). The
#         output will have a start-index array and an end-index array. This way,
#         when traversing the flowline IDs, there will be an array that indicates
#         which indices to use for all upstream flowlines. Index will be 1-based,
#         and provide the start (first) and end (last) index for all upstream flowlines.

#         Flowlines with no upstream segment will have start and end indices of 0.
#         '''

#         newDim = 'index'
#         if addDim:
#             rootgrp.createDimension(newDim, len(FromList))                      # Create Dimension

#         # Build the new from and index variables
#         From = rootgrp.createVariable(fromVar, outDtype, newDim)                # Variable (32-bit signed integer)
#         StartVar = rootgrp.createVariable('fromIndsStart', outDtype, dim)       # Variable (32-bit signed integer)
#         EndVar = rootgrp.createVariable('fromIndsEnd', outDtype, dim)           # Variable (32-bit signed integer)
#         ToVar = rootgrp.createVariable(toIndexVar, outDtype, dim)               # Variable (32-bit signed integer)

#         # Add variable attributes
#         #From.long_name = 'Vector of upstream segments. For use with upstart and upend index variables.'
#         From.long_name = '1-based index into link variable. Represents index of the upstream segment for each flowline in link variables.'
#         StartVar.long_name = '1-based index into %s variable. Represents index of first upstream segment for each flowline in link variable' %(fromVar)
#         EndVar.long_name = '1-based index into %s variable. Represents index of last upstream segment for each flowline in link variable' %(fromVar)
#         ToVar.long_name = '1-based index into to variable. Represents index of the downstream segment for each flowline in link variable'

#         # Iterate over the 'link' vector, finding all upstreams
#         FromArr = numpy.array(FromList)                                         # Array of upstream segments for each downstream segment in FromwhichArr
#         FromwhichArr = numpy.array(FromwhichList)                               # Array of downstream segments for each upstream segment in FromArr

#         # Set up start vector
#         startArr = numpy.full(listlen, 0, dtype=numpy.int)                     # Initialize array with zeros
#         endArr = numpy.full(listlen, 0, dtype=numpy.int)                       # Initialize array with zeros

#         # Find common elements
#         commons = numpy.in1d(idsArr, FromwhichArr)                              # Find common elements between link and from arrays

#         # This block will sort the Fromwhich array and find indices into idsArr
#         index = numpy.argsort(FromwhichArr, kind="mergesort")                   # Find sort order of Fromwhich array. Mergesort preserves order (stable sort)
#         sorted_x = FromwhichArr[index]                                          # Sort the Fromwhich array
#         sorted_index = numpy.searchsorted(sorted_x, idsArr)                     # Find the indices in the sorted Fromwhich array to insert unsorted idsArr
#         yindex = numpy.take(index, sorted_index, mode="clip")                   # Obtain the indices and order to make Fromwhich array match idsArr
#         startArr[commons] = yindex[commons]                                     # Record these indices as the start order
#         del index, sorted_x, sorted_index

#         # Obtain the indices of start and stops of uqnique values in the Fromwhich array
#         # From https://stackoverflow.com/questions/47495510/numpy-in-a-sorted-list-find-the-first-and-the-last-index-for-each-unique-value?rq=1
#         _,ind1,inv1,cou1 = numpy.unique(FromwhichArr, return_index=True, return_inverse=True, return_counts=True)
#         countsArr = cou1[inv1]          # This is the number of unique values for each item along the Fromwhicharr

#         # Diagnostics to see how many tributaries each link has
#         numTribsArr = countsArr[yindex][commons]
#         for val in numpy.unique(numTribsArr):
#             print('  %s links have %s tributaries' %((numTribsArr==val).sum(), val))

#         # Set the end array
#         endArr[commons] = startArr[commons] + numTribsArr-1                     # The -1 will ensure the end index represents the last index to use for from segments
#         del numTribsArr

#         # Block to build index on idsArr for each To segement.
#         tosArr = rootgrp.variables[TOvar][:]                                    # Read the existing 'to' variable from RouteLink
#         index2 = numpy.argsort(idsArr, kind="mergesort")                        # Find the sort order of the idsArr array. Mergesort preserves order (stable sort)
#         sorted_x2 = idsArr[index2]                                              # Sort the idsArr
#         sorted_index2 = numpy.searchsorted(sorted_x2, tosArr)                   # Find the indices in the sorted idsArr array to insert unsorted to-segments
#         yindex2 = numpy.take(index2, sorted_index2, mode="clip")                # Obtain the indices and order to make idsArr match the toArr
#         yindex2[numpy.where(tosArr==0)] = 0                                     # Set all segments with no downstream segment to 0
#         if (tosArr[tosArr!=0] == idsArr[yindex2][tosArr!=0]).sum() == tosArr[tosArr!=0].shape[0]:
#             print('The commons elements between arrays match for to array')
#         #if (idsArr[yindex2]==tosArr).sum() == tosArr.shape:
#         #    print('Using the indices results in recreating the to array')
#         del sorted_index2


#         # Block to build index on idsArr for each From segment in order to populate the From indices variable
#         sorted_index3 = numpy.searchsorted(sorted_x2, FromArr)                   # Find the indices in the sorted idsArr array to insert unsorted to-segments
#         yindex3 = numpy.take(index2, sorted_index3, mode="clip")                # Obtain the indices and order to make idsArr match the toArr
#         if (FromArr == idsArr[yindex3]).sum() == FromArr.shape[0]:
#             # Proves that if you reindex the 'link' variable using FromIndex values, you get the From array
#             print('Using the indices results in recreating the From array')
#         del index2, sorted_x2, sorted_index3

#         # Test to make sure indices match up.
#         if (idsArr[commons] == FromwhichArr[startArr[commons]]).sum() == commons.sum():
#             print('The commons elements between arrays match for start array')
#         if (idsArr[commons] == FromwhichArr[endArr[commons]]).sum() == commons.sum():
#             print('The commons elements between arrays match for end array')

#         # Make the arrays 1-based
#         startArr[commons] += 1                                                  # Make the array 1-based
#         endArr[commons] += 1                                                    # Make the array 1-based
#         yindex2[numpy.where(tosArr!=0)] += 1                                    # Make the array 1-based
#         #yindex2[numpy.where(yindex2>=yindex2.shape[0])] += -1                   # Subtract one if the largest index is out of range for the array size
#         yindex3 += 1

#         # Write to output file
#         #From[:] = FromArr                                                      # This is simply a list of all the from segment ComIDs
#         From[:] = yindex3                                                       # This is a list of all the from indices into the link variable
#         EndVar[:] = endArr
#         StartVar[:] = startArr
#         ToVar[:] = yindex2

#         del tosArr

#     # Write to output variables and close file
#     checkIndices(rootgrp, Tos)                                   # Check file
#     rootgrp.close()
#     del Tos
#     print('        Done writing {0} table to disk.'.format(inNC))
#     print('Finished processing in %3.2f seconds.' %(time.time()-tic))
    
