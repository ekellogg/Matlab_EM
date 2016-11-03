# Pymol script for representing the sectors of 3TGI.

delete all
load /Users/rama/Documents/software/SCA5_forDist/Inputs/3TGI.pdb, main
hide all
bg_color white
show cartoon, (chain E)

color white

create sector_1, (resi 16,18,22,29,54,67,92,94,100,101,108,111,120,123,134,138,141,161,172,176,177,180,183,184,187,188,188A,189,190,192,203,209,213,214,215,216,219,221,222,223,226,227,228,230)& (chain E)

show spheres, sector_1

set_color colres16, [1,0,0]
color colres16, (resi 16)& (chain E)

set_color colres18, [1,0,0]
color colres18, (resi 18)& (chain E)

set_color colres22, [1,0,0]
color colres22, (resi 22)& (chain E)

set_color colres29, [1,0,0]
color colres29, (resi 29)& (chain E)

set_color colres54, [1,0,0]
color colres54, (resi 54)& (chain E)

set_color colres67, [1,0,0]
color colres67, (resi 67)& (chain E)

set_color colres92, [1,0,0]
color colres92, (resi 92)& (chain E)

set_color colres94, [1,0,0]
color colres94, (resi 94)& (chain E)

set_color colres100, [1,0,0]
color colres100, (resi 100)& (chain E)

set_color colres101, [1,0,0]
color colres101, (resi 101)& (chain E)

set_color colres108, [1,0,0]
color colres108, (resi 108)& (chain E)

set_color colres111, [1,0,0]
color colres111, (resi 111)& (chain E)

set_color colres120, [1,0,0]
color colres120, (resi 120)& (chain E)

set_color colres123, [1,0,0]
color colres123, (resi 123)& (chain E)

set_color colres134, [1,0,0]
color colres134, (resi 134)& (chain E)

set_color colres138, [1,0,0]
color colres138, (resi 138)& (chain E)

set_color colres141, [1,0,0]
color colres141, (resi 141)& (chain E)

set_color colres161, [1,0,0]
color colres161, (resi 161)& (chain E)

set_color colres172, [1,0,0]
color colres172, (resi 172)& (chain E)

set_color colres176, [1,0,0]
color colres176, (resi 176)& (chain E)

set_color colres177, [1,0,0]
color colres177, (resi 177)& (chain E)

set_color colres180, [1,0,0]
color colres180, (resi 180)& (chain E)

set_color colres183, [1,0,0]
color colres183, (resi 183)& (chain E)

set_color colres184, [1,0,0]
color colres184, (resi 184)& (chain E)

set_color colres187, [1,0,0]
color colres187, (resi 187)& (chain E)

set_color colres188, [1,0,0]
color colres188, (resi 188)& (chain E)

set_color colres188A, [1,0,0]
color colres188A, (resi 188A)& (chain E)

set_color colres189, [1,0,0]
color colres189, (resi 189)& (chain E)

set_color colres190, [1,0,0]
color colres190, (resi 190)& (chain E)

set_color colres192, [1,0,0]
color colres192, (resi 192)& (chain E)

set_color colres203, [1,0,0]
color colres203, (resi 203)& (chain E)

set_color colres209, [1,0,0]
color colres209, (resi 209)& (chain E)

set_color colres213, [1,0,0]
color colres213, (resi 213)& (chain E)

set_color colres214, [1,0,0]
color colres214, (resi 214)& (chain E)

set_color colres215, [1,0,0]
color colres215, (resi 215)& (chain E)

set_color colres216, [1,0,0]
color colres216, (resi 216)& (chain E)

set_color colres219, [1,0,0]
color colres219, (resi 219)& (chain E)

set_color colres221, [1,0,0]
color colres221, (resi 221)& (chain E)

set_color colres222, [1,0,0]
color colres222, (resi 222)& (chain E)

set_color colres223, [1,0,0]
color colres223, (resi 223)& (chain E)

set_color colres226, [1,0,0]
color colres226, (resi 226)& (chain E)

set_color colres227, [1,0,0]
color colres227, (resi 227)& (chain E)

set_color colres228, [1,0,0]
color colres228, (resi 228)& (chain E)

set_color colres230, [1,0,0]
color colres230, (resi 230)& (chain E)

show surface, sector_1

set transparency, 0.4
create sector_2, (resi 21,24,25,26,29,44,46,52,54,56,66,68,69,71,72,77,81,92,104,105,107,108,118,123,124,134,136,141,153,156,157,180,183,184,199,201,203,210,225,227,229,242,245)& (chain E)

show spheres, sector_2

set_color colres21, [0,0,1]
color colres21, (resi 21)& (chain E)

set_color colres24, [0,0,1]
color colres24, (resi 24)& (chain E)

set_color colres25, [0,0,1]
color colres25, (resi 25)& (chain E)

set_color colres26, [0,0,1]
color colres26, (resi 26)& (chain E)

set_color colres29, [0,0,1]
color colres29, (resi 29)& (chain E)

set_color colres44, [0,0,1]
color colres44, (resi 44)& (chain E)

set_color colres46, [0,0,1]
color colres46, (resi 46)& (chain E)

set_color colres52, [0,0,1]
color colres52, (resi 52)& (chain E)

set_color colres54, [0,0,1]
color colres54, (resi 54)& (chain E)

set_color colres56, [0,0,1]
color colres56, (resi 56)& (chain E)

set_color colres66, [0,0,1]
color colres66, (resi 66)& (chain E)

set_color colres68, [0,0,1]
color colres68, (resi 68)& (chain E)

set_color colres69, [0,0,1]
color colres69, (resi 69)& (chain E)

set_color colres71, [0,0,1]
color colres71, (resi 71)& (chain E)

set_color colres72, [0,0,1]
color colres72, (resi 72)& (chain E)

set_color colres77, [0,0,1]
color colres77, (resi 77)& (chain E)

set_color colres81, [0,0,1]
color colres81, (resi 81)& (chain E)

set_color colres92, [0,0,1]
color colres92, (resi 92)& (chain E)

set_color colres104, [0,0,1]
color colres104, (resi 104)& (chain E)

set_color colres105, [0,0,1]
color colres105, (resi 105)& (chain E)

set_color colres107, [0,0,1]
color colres107, (resi 107)& (chain E)

set_color colres108, [0,0,1]
color colres108, (resi 108)& (chain E)

set_color colres118, [0,0,1]
color colres118, (resi 118)& (chain E)

set_color colres123, [0,0,1]
color colres123, (resi 123)& (chain E)

set_color colres124, [0,0,1]
color colres124, (resi 124)& (chain E)

set_color colres134, [0,0,1]
color colres134, (resi 134)& (chain E)

set_color colres136, [0,0,1]
color colres136, (resi 136)& (chain E)

set_color colres141, [0,0,1]
color colres141, (resi 141)& (chain E)

set_color colres153, [0,0,1]
color colres153, (resi 153)& (chain E)

set_color colres156, [0,0,1]
color colres156, (resi 156)& (chain E)

set_color colres157, [0,0,1]
color colres157, (resi 157)& (chain E)

set_color colres180, [0,0,1]
color colres180, (resi 180)& (chain E)

set_color colres183, [0,0,1]
color colres183, (resi 183)& (chain E)

set_color colres184, [0,0,1]
color colres184, (resi 184)& (chain E)

set_color colres199, [0,0,1]
color colres199, (resi 199)& (chain E)

set_color colres201, [0,0,1]
color colres201, (resi 201)& (chain E)

set_color colres203, [0,0,1]
color colres203, (resi 203)& (chain E)

set_color colres210, [0,0,1]
color colres210, (resi 210)& (chain E)

set_color colres225, [0,0,1]
color colres225, (resi 225)& (chain E)

set_color colres227, [0,0,1]
color colres227, (resi 227)& (chain E)

set_color colres229, [0,0,1]
color colres229, (resi 229)& (chain E)

set_color colres242, [0,0,1]
color colres242, (resi 242)& (chain E)

set_color colres245, [0,0,1]
color colres245, (resi 245)& (chain E)

show surface, sector_2

set transparency, 0.4
create sector_3, (resi 16,17,19,28,30,33,42,43,52,55,56,57,58,66,69,71,95,102,140,141,142,152,155,156,168,180,182,192,193,194,195,196,197,198,199,211,213,214,216,225,227)& (chain E)

show spheres, sector_3

set_color colres16, [0,1,0]
color colres16, (resi 16)& (chain E)

set_color colres17, [0,1,0]
color colres17, (resi 17)& (chain E)

set_color colres19, [0,1,0]
color colres19, (resi 19)& (chain E)

set_color colres28, [0,1,0]
color colres28, (resi 28)& (chain E)

set_color colres30, [0,1,0]
color colres30, (resi 30)& (chain E)

set_color colres33, [0,1,0]
color colres33, (resi 33)& (chain E)

set_color colres42, [0,1,0]
color colres42, (resi 42)& (chain E)

set_color colres43, [0,1,0]
color colres43, (resi 43)& (chain E)

set_color colres52, [0,1,0]
color colres52, (resi 52)& (chain E)

set_color colres55, [0,1,0]
color colres55, (resi 55)& (chain E)

set_color colres56, [0,1,0]
color colres56, (resi 56)& (chain E)

set_color colres57, [0,1,0]
color colres57, (resi 57)& (chain E)

set_color colres58, [0,1,0]
color colres58, (resi 58)& (chain E)

set_color colres66, [0,1,0]
color colres66, (resi 66)& (chain E)

set_color colres69, [0,1,0]
color colres69, (resi 69)& (chain E)

set_color colres71, [0,1,0]
color colres71, (resi 71)& (chain E)

set_color colres95, [0,1,0]
color colres95, (resi 95)& (chain E)

set_color colres102, [0,1,0]
color colres102, (resi 102)& (chain E)

set_color colres140, [0,1,0]
color colres140, (resi 140)& (chain E)

set_color colres141, [0,1,0]
color colres141, (resi 141)& (chain E)

set_color colres142, [0,1,0]
color colres142, (resi 142)& (chain E)

set_color colres152, [0,1,0]
color colres152, (resi 152)& (chain E)

set_color colres155, [0,1,0]
color colres155, (resi 155)& (chain E)

set_color colres156, [0,1,0]
color colres156, (resi 156)& (chain E)

set_color colres168, [0,1,0]
color colres168, (resi 168)& (chain E)

set_color colres180, [0,1,0]
color colres180, (resi 180)& (chain E)

set_color colres182, [0,1,0]
color colres182, (resi 182)& (chain E)

set_color colres192, [0,1,0]
color colres192, (resi 192)& (chain E)

set_color colres193, [0,1,0]
color colres193, (resi 193)& (chain E)

set_color colres194, [0,1,0]
color colres194, (resi 194)& (chain E)

set_color colres195, [0,1,0]
color colres195, (resi 195)& (chain E)

set_color colres196, [0,1,0]
color colres196, (resi 196)& (chain E)

set_color colres197, [0,1,0]
color colres197, (resi 197)& (chain E)

set_color colres198, [0,1,0]
color colres198, (resi 198)& (chain E)

set_color colres199, [0,1,0]
color colres199, (resi 199)& (chain E)

set_color colres211, [0,1,0]
color colres211, (resi 211)& (chain E)

set_color colres213, [0,1,0]
color colres213, (resi 213)& (chain E)

set_color colres214, [0,1,0]
color colres214, (resi 214)& (chain E)

set_color colres216, [0,1,0]
color colres216, (resi 216)& (chain E)

set_color colres225, [0,1,0]
color colres225, (resi 225)& (chain E)

set_color colres227, [0,1,0]
color colres227, (resi 227)& (chain E)

show surface, sector_3

set transparency, 0.4
show stick, (resi 15,16,17) & (chain I)
color yellow, (resi 15,16,17) & (chain I)
set stick_radius, 0.4
show stick, (chain B)
util.cbay main and chain I
color red, sector_1
color blue, sector_2
color green, sector_3
set two_sided_lighting, 1
set sphere_quality, 2
set surface_quality, 2
set stick_quality, 10
set cartoon_oval_quality, 10
set cartoon_loop_quality, 6
set_view (\
-0.655116498,   -0.695889294,   -0.294215441,\
0.386854321,    0.025535638,   -0.921788275,\
0.648974240,   -0.717696548,    0.252478868,\
0.000000000,    0.000000000, -148.182495117,\
-12.619968414,  -89.585556030,   -3.959166050,\
122.692764282,  173.672241211,  -20.000000000 )
