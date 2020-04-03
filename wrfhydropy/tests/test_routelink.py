import numpy as np
import pathlib
import xarray as xr

from wrfhydropy.core.routelink import Routelink

rl_file = pathlib.Path(
    '/Users/james/Downloads/croton_Route_Link_fromVars.nc')
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
    'headwater': [6228190],
    'outlet': [6227150],
    'both': [6228190, 6227150]}

inds_check_answer = {
    'headwater': [20],
    'outlet': [179],
    'both': [20, 179]}

# ---------------------------------------------------------------------------
# id (link) to ind (feature_id) translation
# pass scalars
inds_check = {
    key: [rl.routelink.id_to_ind(intg) for intg in lst]
    for key, lst in ids_check_answer.items()}
for key in inds_check:
    assert all([inds_check[key][ii] == inds_check_answer[key][ii]
                for ii in range(len(inds_check_answer[key]))])

# pass lists
inds_check = {
    key: rl.routelink.ids_to_inds(value)
    for key, value in ids_check_answer.items()}
for key in inds_check:
    assert inds_check[key] == inds_check_answer[key]

# Full domain list
inds_check_answer_full_domain = rl.feature_id.values.tolist()
inds_check_list_full_domain = (
    rl.routelink.ids_to_inds(rl.link.values.tolist()))
assert inds_check_answer_full_domain == inds_check_list_full_domain

# id (link) to id (feature_id) translation
# pass scalars
ids_check = {
    key: [rl.routelink.ind_to_id(intg) for intg in lst]
    for key, lst in inds_check_answer.items()}
for key in ids_check:
    assert all([ids_check[key][ii] == ids_check_answer[key][ii]
                for ii in range(len(ids_check_answer[key]))])

# pass lists
ids_check = {
    key: rl.routelink.inds_to_ids(value)
    for key, value in inds_check_answer.items()}
for key in ids_check:
    assert ids_check[key] == ids_check_answer[key]

# Full domain list
ids_check_answer_full_domain = rl.link.values.tolist()
ids_check_list_full_domain = (
    rl.routelink.inds_to_ids(rl.feature_id.values.tolist()))
assert ids_check_answer_full_domain == ids_check_list_full_domain


# ---------------------------------------------------------------------------
# ind tracing -------------------------------------------------------------------

# Confirm the downstream indices
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
        for key, value in ids_check.items()}
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
        for key, value in ids_check.items()}
    for key in down_ind_id_in_check.keys():
        assert (down_ind_id_in_check[key] ==
                down_ind_check_max_depth_answer[max_depth][key])


# Confirm the upstream indices
# TODO: check birfuractions with max_depth
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


# ---------------------------------------------------------------------------
# id tracing -----------------------------------------------------------
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

for max_depth in down_id_check_max_depth_answer.keys():
    # list ind passed
    down_id_check = {
        key: rl.routelink.get_downstream_ids(value, max_depth=max_depth)
        for key, value in inds_check_answer.items()}
    for key in down_id_check.keys():
        assert (down_id_check[key] ==
                down_id_check_max_depth_answer[max_depth][key])
    # list id passed
    down_id_id_in_check = {
        key: rl.routelink.get_downstream_ids(
            value, max_depth=max_depth, id_in=True)
        for key, value in ids_check.items()}    
    for key in down_id_id_in_check.keys():
        assert (down_id_id_in_check[key] ==
                down_id_check_max_depth_answer[max_depth][key])

# ---------------------------------------------------------------------------
# Get outlet
headwater_ind = inds_check_answer['headwater'][0]
outlet_ind = inds_check_answer['outlet']
headwater_id = ids_check['headwater']
outlet_id = ids_check['outlet']

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

# list input
# Get all outlet gages
all_outlets = rl.routelink.get_outlet_inds(rl.feature_id.values.tolist())
outlet_inds = np.unique([tt[1][0] for tt in all_outlets])
assert (rl.to.isel(feature_id=outlet_inds).values == 0).all()
# also check that these are all the zeros in the to field
assert len(rl.where(rl.to == 0, drop=True).to) == len(outlet_inds)


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

only_gages_from_ids_drop = rl.routelink.ids_to_gages(all_ids, drop_missing=True)
assert only_gages_from_ids_drop[1] == only_gages_check
only_gages_from_ids_no_drop = rl.routelink.ids_to_gages(
    only_gages_from_ids_drop[0], drop_missing=False)
assert only_gages_from_ids_no_drop == only_gages_check


# ---------------------------------------------------------------------------
# TODO test psasing gage to above functions.

# ---------------------------------------------------------------------------
# Get gages below current gage
just_gages = rl.routelink.inds_to_gages(all_inds)
gage_outlets = rl.routelink.get_outlet_inds(just_gages[0])
outlets = np.unique([go[1][0] for go in gage_outlets]).tolist()

outlet_gages = {}
for oo in outlets:
    outlet_gages[oo] = []
    for go in gage_outlets:
        if go[1][0] == oo:
            outlet_gages[oo] += [go[0][0]]

# DOWNstream approach should be faster
outlet_gages_down_gages = {
    key: {vv: rl.routelink.inds_to_gages(
                  rl.routelink.get_downstream_inds(vv)[1])[0]
          for vv in val}
    for key, val in outlet_gages.items()}
outlet_gages_down_gages = {
    k0: {k1: v1 for k1, v1 in v0.items() if v1 != []}
    for k0, v0 in outlet_gages_down_gages.items()}

# Could probably rely on ordering but that's sketch
# which ever value's map produces all the other values in the list, that's the closest
# just keep that one. map becomes 1-1 (but repeated rhs/values)
outlet_gages_down_gage = {}
 k0, v0 in outlet_gages_down_gages.items():
    outlet_gages_down_gage[k0] = {}
     k1, v1 in v0.items():
        if len(v1) == 1:
            outlet_gages_down_gage[k0][k1] = v1
        else:
             vv in v1:
                if vv in v0.keys():
                    if sorted([vv] + v0[vv]) == sorted(v1):
                        outlet_gages_down_gage[k0][k1] = [vv]
                        break

# invert the dictonary to get upstream of each ind: will be one to many 
down_gages = {k0: np.unique([v1 for k1, v1 in v0.items()]).tolist()
              for k0, v0 in outlet_gages_down_gage.items()}
up_gages = {}
for k0, v0 in outlet_gages_down_gage.items():
    up_gages[k0] = {}
    for gg in down_gages[k0]:
        up_gages[k0][gg] = []
        for k1, v1 in outlet_gages_down_gage[k0].items():
            if gg == v1[0]:
                up_gages[k0][gg] += [k1]

affa

# # function to get outlet  each link
# from_ind_check = to_ind_rl_check_answer
# to_ind_rl_check_answer = {
#     key: np.where(rl.link == value)[0].tolist()
#      key, value in to_link_rl_check.items()}
# to_ind_check = {
#     key: rl.routelink.get_downstream_inds(value[0], max_depth=1)[1]
#      key, value in inds_check_scalar_answer.items()}
#  key in to_ind_check.keys():
#     assert np.array_equal(
#         np.array(to_ind_check[key]), np.array(to_ind_rl_check_answer[key]) )


# # -----------------------------------------------------------------------------

# all_downstream_inds = [
#     rl.routelink.get_downstream_inds(link, max_depth=1)  link in upstream_inds[1]]


# len(to_ind_check['outlet'])


# # to re

# trace_down_1 = {
#     key: rl.routelink.get_downstream_inds(value[0], max_depth=1)
#      key, value in inds_to_check.items()}


# # -----------------------------------------------------------------------------

# downstream_inds = rl.routelink.get_downstream_inds(33528)

# # Go back up but just a bit further, bifurcating to the mainstem
# upstream_inds = rl.routelink.get_upstream_inds(downstream_inds[1][0], max_depth=2)
# # Start from the result of the previous, go back down but just one link
# all_downstream_inds = [
#     rl.routelink.get_downstream_inds(link, max_depth=1)  link in upstream_inds[1]]



# #check_from = {key: rl.to[value].values  key, value in inds_to_check.items()}

# check_to

# upstream_inds = rl.routelink.get_upstream_inds(, max_depth=2)
# # Start from the result of the previous, go back down but just one link


# downstream_inds = rl.routelink.get_downstream_inds(33528)

# # Go back up but just a bit further, bifurcating to the mainstem
# all_downstream_inds = [
#     rl.routelink.get_downstream_inds(link, max_depth=1)  link in upstream_inds[1]]



