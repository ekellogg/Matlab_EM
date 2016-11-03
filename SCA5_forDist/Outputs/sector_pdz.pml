# Pymol script for representing the sectors of 1BE9.

delete all
load /Users/rama/Documents/software/SCA5_forDist/Inputs/1BE9.pdb, main
hide all
bg_color white
show cartoon, (chain A)

color white

create sector_1, (resi 322,323,325,327,329,330,336,347,351,353,359,362,363,364,372,375,376,379,386,388)& (chain A)

show spheres, sector_1

set_color colres322, [0,0,1]
color colres322, (resi 322)& (chain A)

set_color colres323, [0,0,1]
color colres323, (resi 323)& (chain A)

set_color colres325, [0,0,1]
color colres325, (resi 325)& (chain A)

set_color colres327, [0,0,1]
color colres327, (resi 327)& (chain A)

set_color colres329, [0,0,1]
color colres329, (resi 329)& (chain A)

set_color colres330, [0,0,1]
color colres330, (resi 330)& (chain A)

set_color colres336, [0,0,1]
color colres336, (resi 336)& (chain A)

set_color colres347, [0,0,1]
color colres347, (resi 347)& (chain A)

set_color colres351, [0,0,1]
color colres351, (resi 351)& (chain A)

set_color colres353, [0,0,1]
color colres353, (resi 353)& (chain A)

set_color colres359, [0,0,1]
color colres359, (resi 359)& (chain A)

set_color colres362, [0,0,1]
color colres362, (resi 362)& (chain A)

set_color colres363, [0,0,1]
color colres363, (resi 363)& (chain A)

set_color colres364, [0,0,1]
color colres364, (resi 364)& (chain A)

set_color colres372, [0,0,1]
color colres372, (resi 372)& (chain A)

set_color colres375, [0,0,1]
color colres375, (resi 375)& (chain A)

set_color colres376, [0,0,1]
color colres376, (resi 376)& (chain A)

set_color colres379, [0,0,1]
color colres379, (resi 379)& (chain A)

set_color colres386, [0,0,1]
color colres386, (resi 386)& (chain A)

set_color colres388, [0,0,1]
color colres388, (resi 388)& (chain A)

show surface, sector_1

set transparency, 0.4
set stick_radius, 0.4
show stick, (chain B)
util.cbay main and chain B
remove main and chain A and resi 300-308
remove main and chain A and resi 394-end
set two_sided_lighting, 1
set sphere_quality, 2
set surface_quality, 2
set stick_quality, 10
set cartoon_oval_quality, 10
set cartoon_loop_quality, 6
set_view (\
-0.880436897,    0.134605691,    0.454657406,\
0.111319803,    0.990739465,   -0.077750385,\
-0.460912317,   -0.017842097,   -0.887267053,\
0.000000000,    0.000000000, -104.157783508,\
35.284568787,   61.732337952,   29.176891327,\
88.666572571,  119.648986816,  -20.000000000 )
