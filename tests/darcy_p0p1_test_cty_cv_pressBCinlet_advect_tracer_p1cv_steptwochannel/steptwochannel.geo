// define a variable used for characteristic length
lc = 10.0;

// define geometric entities

// define points
p_1a = newp; Point(p_1a) = {   0.0,   0.0, 0.0,lc};
p_2a = newp; Point(p_2a) = { 150.0,   0.0, 0.0,lc};
p_3a = newp; Point(p_3a) = { 300.0,   0.0, 0.0,lc};
p_1b = newp; Point(p_1b) = {   0.0, 150.0, 0.0,lc};
p_2b = newp; Point(p_2b) = { 150.0, 150.0, 0.0,lc};
p_3b = newp; Point(p_3b) = { 300.0, 150.0, 0.0,lc};
p_1c = newp; Point(p_1c) = {   0.0, 300.0, 0.0,lc};
p_2c = newp; Point(p_2c) = { 150.0, 300.0, 0.0,lc};
p_3c = newp; Point(p_3c) = { 300.0, 300.0, 0.0,lc};

// define line's
l_12a = newc; Line(l_12a) = {p_1a,p_2a};
l_23a = newc; Line(l_23a) = {p_2a,p_3a};

l_1ab = newc; Line(l_1ab) = {p_1a,p_1b};
l_2ab = newc; Line(l_2ab) = {p_2a,p_2b};
l_3ab = newc; Line(l_3ab) = {p_3a,p_3b};

l_12b = newc; Line(l_12b) = {p_1b,p_2b};
l_23b = newc; Line(l_23b) = {p_2b,p_3b};

l_1bc = newc; Line(l_1bc) = {p_1b,p_1c};
l_2bc = newc; Line(l_2bc) = {p_2b,p_2c};
l_3bc = newc; Line(l_3bc) = {p_3b,p_3c};

l_12c = newc; Line(l_12c) = {p_1c,p_2c};
l_23c = newc; Line(l_23c) = {p_2c,p_3c};

// define line loop
ll_12ab = newll; Line Loop(ll_12ab) = {l_12a,l_2ab,-l_12b,-l_1ab};
ll_23ab = newll; Line Loop(ll_23ab) = {l_23a,l_3ab,-l_23b,-l_2ab};
ll_12bc = newll; Line Loop(ll_12bc) = {l_12b,l_2bc,-l_12c,-l_1bc};
ll_23bc = newll; Line Loop(ll_23bc) = {l_23b,l_3bc,-l_23c,-l_2bc};

// define surface
s_12ab = news; Plane Surface(s_12ab) = {ll_12ab};
s_23ab = news; Plane Surface(s_23ab) = {ll_23ab};
s_12bc = news; Plane Surface(s_12bc) = {ll_12bc};
s_23bc = news; Plane Surface(s_23bc) = {ll_23bc};

// define physical entities

// define physical lines
Physical Line(7) = {l_12a, l_23a};
Physical Line(8) = {l_3ab, l_3bc};
Physical Line(9) = {l_12c, l_23c};
Physical Line(10) = {l_1ab};
Physical Line(100) = {l_1bc};

// define physical surfaces
Physical Surface(11) ={s_12ab};
Physical Surface(12) ={s_23ab};
Physical Surface(13) ={s_12bc};
Physical Surface(14) ={s_23bc};
