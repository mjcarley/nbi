lc = 0.25 ;

r = 1.3 ;

Point(0) = { 0,  0,  0, lc} ;
Point(1) = { r,  0,  0, lc} ;
Point(2) = { 0,  r,  0, lc} ;
Point(3) = { 0,  0,  r, lc} ;
Point(4) = {-r,  0,  0, lc} ;
Point(5) = { 0, -r,  0, lc} ;
Point(6) = { 0,  0, -r, lc} ;



//+
Line(1) = {0, 3};
//+
Line(2) = {0, 2};
//+
Line(3) = {0, 1};
//+
Circle(4) = {1, 0, 6};
//+
Circle(5) = {6, 0, 4};
//+
Circle(6) = {4, 0, 3};
//+
Circle(7) = {3, 0, 1};
//+
Circle(8) = {1, 0, 2};
//+
Circle(9) = {2, 0, 6};
//+
Circle(10) = {2, 0, 4};
//+
Circle(11) = {2, 0, 3};
//+
Circle(12) = {5, 0, 1};
//+
Circle(13) = {5, 0, 6};
//+
Circle(14) = {5, 0, 4};
//+
Circle(15) = {5, 0, 3};
//+
Curve Loop(1) = {3, 8, -2};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {2, 11, -1};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {1, 7, -3};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {7, -12, 15};
//+
Surface(4) = {-4};
//+
Curve Loop(5) = {12, 4, -13};
//+
Surface(5) = {-5};
//+
Curve Loop(6) = {4, -9, -8};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {15, -6, -14};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {13, 5, -14};
//+
Surface(8) = {-8};
//+
Curve Loop(9) = {5, -10, 9};
//+
Surface(9) = {9};
//+
Curve Loop(10) = {11, -6, -10};
//+
Surface(10) = {-10};
