// tube geometry for visualization of acoustic field around rotor and
// nacelle
len = 3 ;
D = 2.5 ;

// mesh resolution setting
lc = 0.1 ;

Point( 1) = {-len/2,    0,    0, lc} ;
Point( 2) = {-len/2,  D/2,    0, lc} ;
Point( 3) = {-len/2,    0,  D/2, lc} ;
Point( 4) = {-len/2, -D/2,    0, lc} ;
Point( 5) = {-len/2,    0, -D/2, lc} ;
	             
Point( 6) = { len/2,    0,    0, lc} ;
Point( 7) = { len/2,  D/2,    0, lc} ;
Point( 8) = { len/2,    0,  D/2, lc} ;
Point( 9) = { len/2, -D/2,    0, lc} ;
Point(10) = { len/2,    0, -D/2, lc} ;

//+
Circle(1) = {7, 6, 8};
//+
Circle(2) = {8, 6, 9};
//+
Circle(3) = {9, 6, 10};
//+
Circle(4) = {10, 6, 7};
//+
Circle(5) = {2, 1, 3};
//+
Circle(6) = {3, 1, 4};
//+
Circle(7) = {4, 1, 5};
//+
Circle(8) = {5, 1, 2};
//+
Line(9) = {7, 2};
//+
Line(10) = {8, 3};
//+
Line(11) = {9, 4};
//+
Line(12) = {10, 5};
//+
Curve Loop(1) = {9, 5, -10, -1};
//+
Surface(1) = {-1};
//+
Curve Loop(2) = {10, 6, -11, -2};
//+
Surface(2) = {-2};
//+
Curve Loop(3) = {11, 7, -12, -3};
//+
Surface(3) = {-3};
//+
Curve Loop(4) = {4, 9, -8, -12};
//+
Surface(4) = {4};
