import numpy as np
import os
import pathlib
import xarray as xr
from .data import collection_data_download

# Not necessary?
# from wrfhydropy.core.routelink import Routelink


# TODO test passing gage to above functions.
# TODO check birfuractions with max_depth
# TODO add routelink_fromVars.nc to the croton test domain.

# Currently "fromVars" is not in the generic domain/files.
# test_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
# The collection_data gets wiped...
# collection_data_download.download()
# rl_file = test_dir / 'data/collection_data/croton_NY/NWM/DOMAIN'
rl_file = pathlib.Path(
#    '/Users/james/Downloads/croton_Route_Link_fromVars.nc')
    '/glade/u/home/jamesmcc/Downloads/croton_Route_Link_fromVars.nc')
rl = xr.open_dataset(rl_file)
# could shorthand: rl_rl = rl.routelink


# These checks were developed in tandem with visualization in
# pywrfhydro/ipynbs/viz/routelink_map.ipynb


# Try to stick to the nomenclature on the LHS
# ind: index, is also the feature_id dimension of routelink as
#      brought in by xarray. Singular   scalars, plural
#      lists/arrays.
# id: id, identifier. in routelink this is 'link'. Singular
#     scalars, plural  lists/arrays. Note that 'id' is a
#     python intrinsic function so I try to avoid using it alone.
# down: downstream, or "to"
# up: upstream, or "from"


# Perform checks on a headwater and an outlet reach
ids_check_answer = {
    'headwater': np.array([6228190]),
    'outlet': np.array([6227150]),
    'both': np.array([6228190, 6227150])}

inds_check_answer = {
    'headwater': np.array([20]),
    'outlet': np.array([179]),
    'both': np.array([20, 179])}


# ind <-> id translation -------------------------------------------------------------------


def test_id_to_ind_scalar():
    # pass scalars
    inds_check = {
        key: [rl.routelink.ids_to_inds(intg) for intg in lst]
        for key, lst in ids_check_answer.items()}
    for key in inds_check:
        assert all([inds_check[key][ii] == inds_check_answer[key][ii]
                    for ii in range(len(inds_check_answer[key]))])


def test_id_to_ind_list():
    # pass lists
    inds_check = {
        key: rl.routelink.ids_to_inds(value)
        for key, value in ids_check_answer.items()}
    for key in inds_check:
        assert (inds_check[key] == inds_check_answer[key]).all()

    # Full domain list
    inds_check_answer_full_domain = rl.feature_id.values.tolist()
    inds_check_list_full_domain = (
        rl.routelink.ids_to_inds(rl.link.values.tolist()))
    assert (inds_check_answer_full_domain ==
            inds_check_list_full_domain).all()


def test_ind_to_id_scalar():
    ids_check = {
        key: [rl.routelink.inds_to_ids(intg) for intg in arr.tolist()]
        for key, arr in inds_check_answer.items()}
    for key in ids_check:
        assert all([ids_check[key][ii] == ids_check_answer[key][ii]
                    for ii in range(len(ids_check_answer[key]))])


def test_ind_to_id_list():
    ids_check = {
        key: rl.routelink.inds_to_ids(value)
        for key, value in inds_check_answer.items()}
    for key in ids_check:
        assert (ids_check[key] == ids_check_answer[key]).all()

    # Full domain list
    ids_check_answer_full_domain = rl.link.values.tolist()
    ids_check_list_full_domain = (
        rl.routelink.inds_to_ids(rl.feature_id.values.tolist()))
    assert (ids_check_answer_full_domain ==
            ids_check_list_full_domain).all()


# ind tracing -------------------------------------------------------------------
down_ind_check_max_depth_answer = {
    1: {
        'headwater': [([20], [59])],
        'outlet': [([179], [])],
        'both': [([20], [59]), ([179], [])]},
    None: {
        'headwater': [
            ([20],
             [59, 93, 108, 117, 130, 136, 142, 152, 158, 160,
              163, 166, 170, 173, 174, 177, 179])],
        'outlet': [
            ([179], [])],
        'both': [
            ([20],
             [59, 93, 108, 117, 130, 136, 142, 152, 158, 160,
              163, 166, 170, 173, 174, 177, 179]),
            ([179], [])]}}


def test_down_ind_trace():
    for max_depth in down_ind_check_max_depth_answer.keys():
        # scalar ind passed
        down_ind_check = {
            key: [rl.routelink.get_downstream_inds(intg, max_depth=max_depth)
                  for intg in value]
            for key, value in inds_check_answer.items()}
        for key in down_ind_check.keys():
            assert (down_ind_check[key] ==
                    down_ind_check_max_depth_answer[max_depth][key])
        # scalar id passed
        down_ind_id_in_check = {
            key: [rl.routelink.get_downstream_inds(
                      intg, max_depth=max_depth, id_in=True)
                  for intg in value]
            for key, value in ids_check_answer.items()}
        for key in down_ind_id_in_check.keys():
            assert (down_ind_id_in_check[key] ==
                    down_ind_check_max_depth_answer[max_depth][key])
        # list ind passed
        down_ind_check = {
            key: rl.routelink.get_downstream_inds(value, max_depth=max_depth)
            for key, value in inds_check_answer.items()}
        for key in down_ind_check.keys():
            assert (down_ind_check[key] ==
                    down_ind_check_max_depth_answer[max_depth][key])
        # list id passed
        down_ind_id_in_check = {
            key: rl.routelink.get_downstream_inds(
                value,max_depth=max_depth,id_in=True)
            for key, value in ids_check_answer.items()}
        for key in down_ind_id_in_check.keys():
            assert (down_ind_id_in_check[key] ==
                    down_ind_check_max_depth_answer[max_depth][key])


# Confirm the upstream indices
up_ind_check_max_depth_answer = {
    1: {'headwater': [([20], [])],
        'outlet': [([179], [177])],
        'both': [([20], []), ([179], [177])]},
    None: {
        'headwater': [([20], [])],
        'outlet': [
            ([179],
             [177, 131, 174, 77, 121, 172, 173, 23, 110, 169,
              109, 170, 89, 167, 87, 166, 78, 155, 164, 76,
              66, 163, 22, 149, 96, 161, 21, 34, 64, 160,
              146, 67, 157, 30, 158, 139, 35, 154, 152, 153,
              128, 148, 142, 150, 151, 122, 62, 143, 126, 136,
              144, 145, 111, 29, 55, 135, 120, 130, 138, 140,
              99, 24, 129, 74, 106, 117, 127, 133, 17, 58,
              73, 119, 14, 86, 108, 61, 118, 116, 15, 16,
              13, 104, 54, 75, 93, 28, 107, 51, 101, 11,
              84, 9, 19, 59, 91, 92, 43, 98, 80, 18,
              20, 56, 57, 48, 4, 25, 27, 44])],
        'both': [
            ([20], []),
            ([179],
             [177, 131, 174, 77, 121, 172, 173, 23, 110, 169,
              109, 170, 89, 167, 87, 166, 78, 155, 164, 76,
              66, 163, 22, 149, 96, 161, 21, 34, 64, 160,
              146, 67, 157, 30, 158, 139, 35, 154, 152, 153,
              128, 148, 142, 150, 151, 122, 62, 143, 126, 136,
              144, 145, 111, 29, 55, 135, 120, 130, 138, 140,
              99, 24, 129, 74, 106, 117, 127, 133, 17, 58,
              73, 119, 14, 86, 108, 61, 118, 116, 15, 16,
              13, 104, 54, 75, 93, 28, 107, 51, 101, 11,
              84, 9, 19, 59, 91, 92, 43, 98, 80, 18,
              20, 56, 57, 48, 4, 25, 27, 44])]} }

def test_up_ind_trace():
    for max_depth in up_ind_check_max_depth_answer.keys():
        up_ind_check = {
            key: rl.routelink.get_upstream_inds(
                value, max_depth=max_depth)
            for key, value in inds_check_answer.items()}
        for key in up_ind_check.keys():
            assert (up_ind_check[key] ==
                    up_ind_check_max_depth_answer[max_depth][key])
        up_ind_id_in_check = {
            key: rl.routelink.get_upstream_inds(
                value, max_depth=max_depth, id_in=True)
            for key, value in ids_check_answer.items()}
        for key in up_ind_id_in_check.keys():
            assert (up_ind_id_in_check[key] ==
                    up_ind_check_max_depth_answer[max_depth][key])


# id tracing -----------------------------------------------------------
# Only testing this downstream since it really not dependent on direction,
# mostly testing the interface, teh guts used the id <-> ind translation
# tested above.

# Translate the previous test from ind to id
down_id_check_max_depth_answer = {}
for key_1, val_1 in down_ind_check_max_depth_answer.items():
    down_id_check_max_depth_answer[key_1] = {}
    for key_2, val_2 in val_1.items():
        down_id_check_max_depth_answer[key_1][key_2] = []
        for ii in range(len(val_2)):
            down_id_check_max_depth_answer[key_1][key_2] += [
                (rl.routelink.inds_to_ids(val_2[ii][0]),
                 rl.routelink.inds_to_ids(val_2[ii][1]))]


def check_down_id(to_check, max_depth):
    answer = down_id_check_max_depth_answer
    for key in to_check.keys():
        for pp in range(len(list(to_check[key]))):
            assert (to_check[key][pp][0] ==
                    answer[max_depth][key][pp][0]).all()
            if len(to_check[key][pp][1]) is 0:
                assert len(answer[max_depth][key][pp][1]) is 0
            else:
                assert (to_check[key][pp][1] ==
                        answer[max_depth][key][pp][1]).all()
    return None


def test_down_id_trace():
    for max_depth in down_id_check_max_depth_answer.keys():
        # list ind passed
        down_id_check = {
            key: rl.routelink.get_downstream_ids(value, max_depth=max_depth)
            for key, value in inds_check_answer.items()}
        _ = check_down_id(down_id_check, max_depth)
        down_id_id_in_check = {
            key: rl.routelink.get_downstream_ids(
                value, max_depth=max_depth, id_in=True)
            for key, value in ids_check_answer.items()}
        _ = check_down_id(down_id_id_in_check, max_depth)


# ---------------------------------------------------------------------------
# Get outlet
headwater_ind = inds_check_answer['headwater'][0]
outlet_ind = inds_check_answer['outlet']
headwater_id = ids_check_answer['headwater']
outlet_id = ids_check_answer['outlet']


def test_get_outlet_inds_ids_scalar():
    # scalar inputs
    # headwater to outlet
    assert (rl.routelink.get_outlet_inds(headwater_ind) ==
            [([headwater_ind], outlet_ind)])
    assert (rl.routelink.get_outlet_inds(headwater_id, id_in=True) ==
            [([headwater_ind], outlet_ind)])
    assert (rl.routelink.get_outlet_ids(headwater_ind) ==
            [(headwater_id, outlet_id)])
    assert (rl.routelink.get_outlet_ids(headwater_id, id_in=True) ==
            [(headwater_id, outlet_id)])
    # outlet to itself
    assert (rl.routelink.get_outlet_inds(outlet_ind) ==
            [(outlet_ind, outlet_ind)])
    assert (rl.routelink.get_outlet_inds(outlet_id, id_in=True) ==
            [(outlet_ind, outlet_ind)])
    assert (rl.routelink.get_outlet_ids(outlet_ind) ==
            [(outlet_id, outlet_id)])
    assert (rl.routelink.get_outlet_ids(outlet_id, id_in=True) ==
            [(outlet_id, outlet_id)])


def test_get_outlet_inds_ids_list():
    # list input
    # Get all outlet gages
    all_outlets = rl.routelink.get_outlet_inds(rl.feature_id.values.tolist())
    outlet_inds = np.unique([tt[1][0] for tt in all_outlets])
    assert (rl.to.isel(feature_id=outlet_inds).values == 0).all()
    # also check that these are all the zeros in the to field
    assert len(rl.where(rl.to == 0, drop=True).to) == len(outlet_inds)


# ---------------------------------------------------------------------------
# Get gages from inds
def test_gages_from_inds():
    all_inds = rl.feature_id.values.tolist()
    all_ids = rl.link.values.tolist()
    all_gages_check = rl.gages.values.tolist()
    all_gages_from_inds = rl.routelink.inds_to_gages(all_inds, drop_missing=False)
    assert all_gages_from_inds == all_gages_check
    all_gages_from_ids = rl.routelink.ids_to_gages(all_ids, drop_missing=False)
    assert all_gages_from_ids == all_gages_check

    missing_gage = b'               '
    only_gages_check = rl.gages.where(rl.gages != missing_gage, drop=True).values.tolist()
    only_gages_from_inds_drop = rl.routelink.inds_to_gages(all_inds, drop_missing=True)
    assert only_gages_from_inds_drop[1] == only_gages_check
    only_gages_from_inds_no_drop = rl.routelink.inds_to_gages(
        only_gages_from_inds_drop[0], drop_missing=False)
    assert only_gages_from_inds_no_drop == only_gages_check

    only_gages_from_ids_drop = rl.routelink.ids_to_gages(all_ids, drop_missing=True)
    assert only_gages_from_ids_drop[1] == only_gages_check
    only_gages_from_ids_no_drop = rl.routelink.ids_to_gages(
        only_gages_from_ids_drop[0], drop_missing=False)
    assert only_gages_from_ids_no_drop == only_gages_check


def test_inds_ids_from_gages():
    all_gages = rl.gages.values.tolist()
    all_inds = rl.feature_id.values.tolist()
    all_ids = rl.link.values.tolist()
    #                123456789012345
    missing_gage = b'               '    
    check_ids = [ii for ii, gg in zip(all_ids, all_gages) if gg != missing_gage]
    check_inds = [ii for ii, gg in zip(all_inds, all_gages) if gg != missing_gage]    
    all_ids_from_gages = rl.routelink.gages_to_ids(all_gages)
    all_inds_from_gages = rl.routelink.gages_to_inds(all_gages)
    for ii in range(len(check_inds)):
        assert all_inds_from_gages[ii] == check_inds[ii]
        assert all_ids_from_gages[check_inds[ii]] == check_ids[ii]

# ---------------------------------------------------------------------------
# Get gages from inds
nested_gage_inds_answer = {179: {160: [142], 179: [143, 160]}}
nested_gages_answer = {
    179: {160: [b'       01374559'],
          179: [b'       01374598', b'       01374581']}}


def test_nested_gage_inds():
    nested_gage_inds = rl.routelink.get_nested_gage_inds()
    assert nested_gage_inds == nested_gage_inds_answer


def test_nested_gages():
    nested_gages = rl.routelink.get_nested_gages()
    assert nested_gages == nested_gages_answer
