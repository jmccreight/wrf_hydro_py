import numpy as np
import pathlib
import xarray as xr

from wrfhydropy.core.routelink import Routelink

rl_file = pathlib.Path(
    '/Users/james/Downloads/croton_Route_Link_fromVars.nc')
rl = xr.open_dataset(rl_file)

# These checks were developed in tandem with visualization in
# pywrfhydro/ipynbs/viz/routelink_map.ipynb

# Try to stick to the nomenclature on the LHS
# ind: index, is also the feature_id dimension of routelink as
#      brought in by xarray. Singular for  scalars, plural for
#      lists/arrays.
# id: id, identifier. in routelink this is 'link'. Singular for
#     scalars, plural for lists/arrays. Note that 'id' is a
#     python intrinsic function so I try to avoid using it alone.
# down: downstream, or "to"
# up: upstream, or "from"

# Perform checks on a headwater and an outlet reach
ids_check = {
    'headwater': 6228190,
    'outlet': 6227150}

# ---------------------------------------------------------------------------
# id (link) to ind (feature_id) translation
inds_check_scalar_answer = {
    'headwater': 20,
    'outlet': 179}
inds_check_scalar = {
    key: rl.routelink.id_to_ind(value)
    for key, value in ids_check.items()}
for key in inds_check_scalar:
    assert inds_check_scalar[key] == inds_check_scalar_answer[key]

inds_check_list_answer = rl.feature_id.values.tolist()
inds_check_list = rl.routelink.ids_to_inds(rl.link.values.tolist())
assert inds_check_list_answer == inds_check_list

# id (link) to ind (feature_id) translation
ids_check_scalar = {
    key: rl.routelink.ind_to_id(value)
    for key, value in inds_check_scalar_answer.items()}
for key in ids_check_scalar:
    assert ids_check_scalar[key] == ids_check[key]

ids_check_list_answer = rl.link.values.tolist()
ids_check_list = rl.routelink.inds_to_ids(rl.feature_id.values.tolist())
assert ids_check_list_answer == ids_check_list


# ---------------------------------------------------------------------------
# ind tracing -------------------------------------------------------------------

# Confirm the downstream indices
down_ind_check_max_depth_answer = {
    1: {'headwater': [59], 'outlet': []},
    None: {'headwater': [
               59, 93, 108, 117, 130, 136, 142, 152, 158, 160,
               163, 166, 170, 173, 174, 177, 179],
           'outlet': []}}
for max_depth in down_ind_check_max_depth_answer.keys():
    down_ind_check = {
        key: rl.routelink.get_downstream_inds(
            value, max_depth=max_depth)[1]  # can this index be removed?
        for key, value in inds_check_scalar_answer.items() }
    for key in down_ind_check.keys():
        assert np.array_equal(
            np.array(down_ind_check[key]),
            np.array(down_ind_check_max_depth_answer[max_depth][key]))
    down_ind_id_in_check = {
        key: rl.routelink.get_downstream_inds(
            value, max_depth=max_depth, id_in=True)[1]  # can this index be removed?
        for key, value in ids_check.items()}
    for key in down_ind_id_in_check.keys():
        assert np.array_equal(
            np.array(down_ind_id_in_check[key]),
            np.array(down_ind_check_max_depth_answer[max_depth][key]))

# Confirm the upstream indices
up_ind_check_max_depth_answer = {
    1: {'headwater': [], 'outlet': [177]},
    None: {'headwater': [],
           'outlet': [
               177, 131, 174, 77, 121, 172, 173, 23, 110, 169,
               109, 170, 89, 167, 87, 166, 78, 155, 164, 76,
               66, 163, 22, 149, 96, 161, 21, 34, 64, 160,
               146, 67, 157, 30, 158, 139, 35, 154, 152, 153,
               128, 148, 142, 150, 151, 122, 62, 143, 126, 136,
               144, 145, 111, 29, 55, 135, 120, 130, 138, 140,
               99, 24, 129, 74, 106, 117, 127, 133, 17, 58,
               73, 119, 14, 86, 108, 61, 118, 116, 15, 16,
               13, 104, 54, 75, 93, 28, 107, 51, 101, 11,
               84, 9, 19, 59, 91, 92, 43, 98, 80, 18,
               20, 56, 57, 48, 4, 25, 27, 44]}}

for max_depth in up_ind_check_max_depth_answer.keys():
    up_ind_check = {
        key: rl.routelink.get_upstream_inds(
            value, max_depth=max_depth)[1]
        for key, value in inds_check_scalar_answer.items()}
    for key in up_ind_check.keys():
        assert np.array_equal(
            np.array(up_ind_check[key]),
            np.array(up_ind_check_max_depth_answer[max_depth][key]))
    up_ind_id_in_check = {
        key: rl.routelink.get_upstream_inds(
            value, max_depth=max_depth, id_in=True)[1]
        for key, value in ids_check.items()}
    for key in up_ind_id_in_check.keys():
        assert np.array_equal(
            np.array(up_ind_id_in_check[key]),
            np.array(up_ind_check_max_depth_answer[max_depth][key]))


# ---------------------------------------------------------------------------
# id tracing -----------------------------------------------------------
# Translate the previous test from ind to id
down_id_check_max_depth_answer = {}
for key_1, val_1 in down_ind_check_max_depth_answer.items():
    down_id_check_max_depth_answer[key_1] = {}
    for key_2, val_2 in val_1.items():
        down_id_check_max_depth_answer[key_1][key_2] = (
            rl.routelink.inds_to_ids(val_2))

for max_depth in down_id_check_max_depth_answer.keys():
    down_id_check = {
        key: rl.routelink.get_downstream_ids(
            value, max_depth=max_depth)
        for key, value in inds_check_scalar.items()}  # this is inds
    for key in down_ind_check.keys():
        assert np.array_equal(
            np.array(down_id_check[key]),
            np.array(down_id_check_max_depth_answer[max_depth][key]))
    down_id_id_in_check = {
        key: rl.routelink.get_downstream_ids(
            value, max_depth=max_depth, id_in=True)
        for key, value in ids_check_scalar.items()}  # this is ids
    for key in down_id_id_in_check.keys():
        assert np.array_equal(
            np.array(down_id_id_in_check[key]),
            np.array(down_id_check_max_depth_answer[max_depth][key]))


# ---------------------------------------------------------------------------
# Get outlet
headwater_ind = inds_check_scalar_answer['headwater']
outlet_ind = inds_check_scalar_answer['outlet']
headwater_id = ids_check['headwater']
outlet_id = ids_check['outlet']
assert rl.routelink.get_outlet_ind(headwater_ind) == outlet_ind
assert rl.routelink.get_outlet_id(headwater_ind) == outlet_id
assert rl.routelink.get_outlet_ind(outlet_ind) == outlet_ind
assert rl.routelink.get_outlet_id(outlet_ind) == outlet_id
assert rl.routelink.get_outlet_ind(headwater_id, id_in=True) == outlet_ind
assert rl.routelink.get_outlet_id(headwater_id, id_in=True) == outlet_id
assert rl.routelink.get_outlet_ind(outlet_id, id_in=True) == outlet_ind
assert rl.routelink.get_outlet_id(outlet_id, id_in=True) == outlet_id

# ---------------------------------------------------------------------------
# Get gages from inds
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

only_gages_from_inds_drop = rl.routelink.ids_to_gages(all_ids, drop_missing=True)
assert only_gages_from_ids_drop[1] == only_gages_check
only_gages_from_ids_no_drop = rl.routelink.ids_to_gages(
    only_gages_from_ids_drop[0], drop_missing=False)
assert only_gages_from_ids_no_drop == only_gages_check



# function to get outlet for each link
from_ind_check = to_ind_rl_check_answer

to_ind_rl_check_answer = {
    key: np.where(rl.link == value)[0].tolist()
    for key, value in to_link_rl_check.items()}
to_ind_check = {
    key: rl.routelink.get_downstream_inds(value[0], max_depth=1)[1]
    for key, value in inds_check_scalar_answer.items()}
for key in to_ind_check.keys():
    assert np.array_equal(
        np.array(to_ind_check[key]), np.array(to_ind_rl_check_answer[key]) )

    
    
# -----------------------------------------------------------------------------

all_downstream_inds = [
    rl.routelink.get_downstream_inds(link, max_depth=1) for link in upstream_inds[1]]


len(to_ind_check['outlet'])


# to re

trace_down_1 = {
    key: rl.routelink.get_downstream_inds(value[0], max_depth=1)
    for key, value in inds_to_check.items()}


# -----------------------------------------------------------------------------

downstream_inds = rl.routelink.get_downstream_inds(33528)

# Go back up but just a bit further, bifurcating to the mainstem
upstream_inds = rl.routelink.get_upstream_inds(downstream_inds[1][0], max_depth=2)
# Start from the result of the previous, go back down but just one link
all_downstream_inds = [
    rl.routelink.get_downstream_inds(link, max_depth=1) for link in upstream_inds[1]]



#check_from = {key: rl.to[value].values for key, value in inds_to_check.items()}

check_to

upstream_inds = rl.routelink.get_upstream_inds(, max_depth=2)
# Start from the result of the previous, go back down but just one link


downstream_inds = rl.routelink.get_downstream_inds(33528)

# Go back up but just a bit further, bifurcating to the mainstem
all_downstream_inds = [
    rl.routelink.get_downstream_inds(link, max_depth=1) for link in upstream_inds[1]]



