rafx = 16;
rafy = 16;
rafz = 1;
Point(1) = {-1, -1, 0, raf};
Point(2) = {1, -1, 0, raf};
Point(3) = {1, 1, 0, raf};
Point(4) = {-1, 1, 0, raf};
Point(5) = {-1, -1, 1, raf};
Point(6) = {1, -1, 1, raf};
Point(7) = {1, 1, 1, raf};
Point(8) = {-1, 1, 1, raf};
Line(1) = {1, 2};
Transfinite Line {1} = rafx+1 Using Progression 1;
Line(2) = {2, 3};
Transfinite Line {2} = rafy+1 Using Progression 1;
Line(3) = {3, 4};
Transfinite Line {3} = rafx+1 Using Progression 1;
Line(4) = {4, 1};
Transfinite Line {4} = rafy+1 Using Progression 1;
Line(5) = {5, 6};
Transfinite Line {5} = rafx+1 Using Progression 1;
Line(6) = {6, 7};
Transfinite Line {6} = rafy+1 Using Progression 1;
Line(7) = {7, 8};
Transfinite Line {7} = rafx+1 Using Progression 1;
Line(8) = {8, 5};
Transfinite Line {8} = rafy+1 Using Progression 1;
Line(9) = {2, 6};
Transfinite Line {9} = rafz+1 Using Progression 1;
Line(10) = {3, 7};
Transfinite Line {10} = rafz+1 Using Progression 1;
Line(11) = {4, 8};
Transfinite Line {11} = rafz+1 Using Progression 1;
Line(12) = {1, 5};
Transfinite Line {12} = rafz+1 Using Progression 1;
Line Loop(1) = {1, 9, -5, -12};
Ruled Surface(1) = {1};
Transfinite Surface {1};
Line Loop(2) = {2, 10, -6, -9};
Ruled Surface(2) = {2};
Transfinite Surface {2};
Line Loop(3) = {3,11,-7,-10};
Ruled Surface(3) = {3};
Transfinite Surface {3};
Line Loop(4) = {4,12,-8,-11};
Ruled Surface(4) = {4};
Transfinite Surface {4};
Line Loop(5) = {-1,-4,-3,-2};
Ruled Surface(5) = {5};
Transfinite Surface {5};
Line Loop(6) = {5,6,7,8};
Ruled Surface(6) = {6};
Transfinite Surface {6};
Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};
Recombine Surface {5};
Recombine Surface {6};
Physical Surface(20) = {1,2,3,4,5,6};
Surface Loop(1) = {6, 1, 5, 4, 3, 2};
Volume(14) = {1};
Transfinite Volume {14};
Physical Volume(1000) = {14};
Mesh.ElementOrder=2 ;
Mesh.SecondOrderIncomplete=1;
