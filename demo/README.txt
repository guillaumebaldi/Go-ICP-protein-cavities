# Coordinates are in PCD format, after < DATA ascii >
# x, y, z, rgb
# Angstrom is the coordinates unit
# Les points ont été déterminés à partir de voxel de taille 1.5 A de coté


# a pair of a same object, one being translated by 10 Angstrom (coordinates unit) along x axis:
2rh1_cavity.pcd		shiftx_2rh1_cavity.pcd

# EASY1 dataset
# pairs of similar cavities (higher registration scores expected):
2c6t_cavity.pcd 	1dm2_cavity.pcd
1ftl_cavity.pcd  	1lb9_cavity.pcd
2b7z_cavity.pcd  	1c6x_cavity.pcd
2ouz_cavity.pcd  	3ert_cavity.pcd
2rh1_cavity.pcd  	5d6l_cavity.pcd

# pairs of dissimilar cavities (lower registration scores expected):
2rh1_cavity.pcd		1dm2_cavity.pcd
2rh1_cavity.pcd		2ouz_cavity.pcd
2rh1_cavity.pcd		2b7z_cavity.pcd
1dm2_cavity.pcd		2ouz_cavity.pcd
1dm2_cavity.pcd		2b7z_cavity.pcd

2c6t_cavity.pcd     145
1dm2_cavity.pcd     129
1ftl_cavity.pcd  	165
1lb9_cavity.pcd     172
2b7z_cavity.pcd  	252
1c6x_cavity.pcd     199
2ouz_cavity.pcd  	217
3ert_cavity.pcd     214
2rh1_cavity.pcd  	170
5d6l_cavity.pcd     158
