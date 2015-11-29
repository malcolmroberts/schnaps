cl__1 = 1e+22;
Point(1) = {-0.7, -0.7, 0, 1e+22};
Point(2) = {0.7, -0.7, 0, 1e+22};
Point(3) = {0.7, 0.7, 0, 1e+22};
Point(4) = {-0.7, 0.7, 0, 1e+22};
Point(5) = {0, 0, 0, 1e+22};
Point(6) = {1.414, 1.414, 0, 1e+22};
Point(7) = {-1.414, 1.414, 0, 1e+22};
Point(8) = {-1.414, -1.414, 0, 1e+22};
Point(9) = {1.414, -1.414, 0, 1e+22};
Point(10) = {0.7, 0.7, 1, 1e+22};
Point(11) = {1.414, 1.414, 1, 1e+22};
Point(15) = {0, 0, 1, 1e+22};
Point(16) = {-1.414, 1.414, 1, 1e+22};
Point(20) = {-0.7, 0.7, 1, 1e+22};
Point(27) = {-1.414, -1.414, 1, 1e+22};
Point(31) = {-0.7, -0.7, 1, 1e+22};
Point(33) = {0.7, -0.7, 1, 1e+22};
Point(35) = {1.414, -1.414, 1, 1e+22};
Line(1) = {1, 2};
Transfinite Line {1} = 9Using Progression 1;
Line(2) = {2, 3};
Transfinite Line {2} = 9Using Progression 1;
Line(3) = {3, 4};
Transfinite Line {3} = 9Using Progression 1;
Line(4) = {4, 1};
Transfinite Line {4} = 9Using Progression 1;
Circle(10) = {9, 5, 6};
Transfinite Line {10} = 9Using Progression 1;
Circle(11) = {6, 5, 7};
Transfinite Line {11} = 9Using Progression 1;
Circle(12) = {7, 5, 8};
Transfinite Line {12} = 9Using Progression 1;
Circle(13) = {8, 5, 9};
Transfinite Line {13} = 9Using Progression 1;
Line(14) = {2, 9};
Transfinite Line {14} = 9Using Progression 1;
Line(15) = {3, 6};
Transfinite Line {15} = 9Using Progression 1;
Line(16) = {4, 7};
Transfinite Line {16} = 9Using Progression 1;
Line(17) = {1, 8};
Transfinite Line {17} = 9Using Progression 1;
Line(27) = {10, 11};
Transfinite Line {27} = 9Using Progression 1;
Circle(28) = {11, 15, 16};
Transfinite Line {28} = 9Using Progression 1;
Line(29) = {16, 20};
Transfinite Line {29} = 9Using Progression 1;
Line(30) = {20, 10};
Transfinite Line {30} = 9Using Progression 1;
Line(32) = {3, 10};
Transfinite Line {32} = 2Using Progression 1;
Line(33) = {6, 11};
Transfinite Line {33} = 2Using Progression 1;
Line(37) = {7, 16};
Transfinite Line {37} = 2Using Progression 1;
Line(41) = {4, 20};
Transfinite Line {41} = 2Using Progression 1;
Circle(50) = {16, 15, 27};
Transfinite Line {50} = 9Using Progression 1;
Line(51) = {27, 31};
Transfinite Line {51} = 9Using Progression 1;
Line(52) = {31, 20};
Transfinite Line {52} = 9Using Progression 1;
Line(59) = {8, 27};
Transfinite Line {59} = 2Using Progression 1;
Line(63) = {1, 31};
Transfinite Line {63} = 2Using Progression 1;
Line(71) = {31, 33};
Transfinite Line {71} = 9Using Progression 1;
Line(72) = {33, 10};
Transfinite Line {72} = 9Using Progression 1;
Line(77) = {2, 33};
Transfinite Line {77} = 2Using Progression 1;
Line(93) = {33, 35};
Transfinite Line {93} = 9Using Progression 1;
Circle(94) = {35, 15, 11};
Transfinite Line {94} = 9Using Progression 1;
Line(99) = {9, 35};
Transfinite Line {99} = 2Using Progression 1;
Circle(116) = {27, 15, 35};
Transfinite Line {116} = 9Using Progression 1;
Line Loop(19) = {14, 10, -15, -2};
Plane Surface(19) = {19};
Transfinite Surface {19};
Recombine Surface {19};
Line Loop(21) = {15, 11, -16, -3};
Plane Surface(21) = {21};
Transfinite Surface {21};
Recombine Surface {21};
Line Loop(23) = {16, 12, -17, -4};
Plane Surface(23) = {23};
Transfinite Surface {23};
Recombine Surface {23};
Line Loop(25) = {17, 13, -14, -1};
Plane Surface(25) = {25};
Transfinite Surface {25};
Recombine Surface {25};
Line Loop(34) = {15, 33, -27, -32};
Ruled Surface(34) = {34};
Transfinite Surface {34};
Recombine Surface {34};
Line Loop(38) = {11, 37, -28, -33};
Ruled Surface(38) = {38};
Transfinite Surface {38};
Recombine Surface {38};
Line Loop(42) = {-16, 41, -29, -37};
Ruled Surface(42) = {42};
Transfinite Surface {42};
Recombine Surface {42};
Line Loop(46) = {-3, 32, -30, -41};
Ruled Surface(46) = {46};
Transfinite Surface {46};
Recombine Surface {46};
Line Loop(47) = {27, 28, 29, 30};
Plane Surface(47) = {47};
Transfinite Surface {47};
Recombine Surface {47};
Line Loop(60) = {12, 59, -50, -37};
Ruled Surface(60) = {60};
Transfinite Surface {60};
Recombine Surface {60};
Line Loop(64) = {-17, 63, -51, -59};
Ruled Surface(64) = {64};
Transfinite Surface {64};
Recombine Surface {64};
Line Loop(68) = {-4, 41, -52, -63};
Ruled Surface(68) = {68};
Transfinite Surface {68};
Recombine Surface {68};
Line Loop(69) = {-29, 50, 51, 52};
Plane Surface(69) = {69};
Transfinite Surface {69};
Recombine Surface {69};
Line Loop(78) = {1, 77, -71, -63};
Ruled Surface(78) = {78};
Transfinite Surface {78};
Recombine Surface {78};
Line Loop(82) = {2, 32, -72, -77};
Ruled Surface(82) = {82};
Transfinite Surface {82};
Recombine Surface {82};
Line Loop(100) = {14, 99, -93, -77};
Ruled Surface(100) = {100};
Transfinite Surface {100};
Recombine Surface {100};
Line Loop(104) = {10, 33, -94, -99};
Ruled Surface(104) = {104};
Transfinite Surface {104};
Recombine Surface {104};
Line Loop(113) = {93, 94, -27, -72};
Plane Surface(113) = {113};
Transfinite Surface {113};
Recombine Surface {113};
Line Loop(126) = {13, 99, -116, -59};
Ruled Surface(126) = {126};
Transfinite Surface {126};
Recombine Surface {126};
Line Loop(135) = {-51, 116, -93, -71};
Plane Surface(135) = {135};
Transfinite Surface {135};
Recombine Surface {135};
Surface Loop(1) = {21, 47, 34, 38, 42, 46};
Volume(1) = {1};
Transfinite Volume {1};
Surface Loop(2) = {23, 69, 42, 60, 64, 68};
Volume(2) = {2};
Transfinite Volume {2};
Surface Loop(4) = {19, 113, 100, 104, 34, 82};
Volume(4) = {4};
Transfinite Volume {4};
Surface Loop(5) = {25, 135, 64, 126, 100, 78};
Volume(5) = {5};
Transfinite Volume {5};
Mesh.ElementOrder=2 ;
Mesh.SecondOrderIncomplete = 1 ;