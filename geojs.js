/**
 * GeoJS
 * (c) 2013 Baimuratov Iskander aka MrBizor
 * 
 * coordinate transformation
 */

 var PI = Math.PI;

function pow(a, c){
	var b = Math.pow(a, c);
	return b;
}

function sqrt(a){
	var b = Math.sqrt(a);
	return b;
}

function toRadians(a){
	var r = a * Math.PI / 180;
	return r;
}

function toDegrees(r){
	var d = r * 180 / Math.PI;
	return d;
}

function sin(a){
	var b = Math.sin(a);
	return b;
}

function cos(a){
	var b = Math.cos(a);
	return b;
}

function tan(a){
	var b = Math.tan(a);
	return b;
}

function atan(a){
	var b = Math.atan(a);
	return b;
}

function trunc(a){
	var b = Math.floor(a);
	return b;
}

function sec(a){
	var b = 1 / cos(a);
	return b;
}

function PrefInt(number, len) {
   return (Array(len).join('0') + number).slice(-length);
}

function dmstodeg(dms){
	var sdms = dms.split(' ');
	var d = parseFloat(sdms[0]);
	var m = parseFloat(sdms[1]);
	var s = parseFloat(sdms[2]);
	var deg = d + m / 60 + s /3600;
	return deg;
}

function degtodms(deg){
	var d = trunc(deg);
	var m = trunc((deg - d) * 60);
	var s = (((deg - d) * 60) - m) * 60;
	var dms = d + '&deg ' + PrefInt(m, 2) + "' " + s.toFixed(5)+'" ';
	return dms;
}

////////////////////////////////////////

function datum(datum_set){
	this.name = datum_set.name;
	this.a = datum_set.a;
	this.f = 1 / datum_set.f;
	this.dx = datum_set.dx;
	this.dy = datum_set.dy;
	this.dz = datum_set.dz;
	this.wx = datum_set.wx;
	this.wy = datum_set.wy;
	this.wz = datum_set.wz;
	this.scale = datum_set.scale;
	
	this.b = datum_set.a * (1 - this.f);
	this.e2 = (pow(datum_set.a, 2) - pow(this.b, 2)) / pow(datum_set.a, 2);
	this.e12 = this.e2 / (1 - this.e2);
	this.n = (datum_set.a - this.b) / (datum_set.a + this.b);
	this.G = datum_set.a * (1 - this.n) * (1 - pow(this.n, 2)) * (1 + (9 * pow(this.n, 2) / 4) +(225 * pow(this.n, 4) / 64)) * (PI / 180);
	this.G0 = 1 - (this.e2/4) - (3 * pow(this.e2, 2) / 64) - (5 * pow(this.e2, 3) / 256);
	this.G2 = 3.0 / 8.0 * (this.e2 + pow(this.e2, 2) /4 + 15 * pow(this.e2, 3) / 128);
	this.G4 = 15.0 / 256.0 * (pow(this.e2, 2) + 3 * pow(this.e2, 3) / 4);
	this.G6 = 35 * pow(this.e2, 3) / 3072;
}

datum.prototype.merdist = function(B){
	var md = this.a * (this.G0 * B - this.G2 * sin (2 * B) + this.G4 * sin (4 * B) - this.G6 * sin (6 * B));
	return md;
}

datum.prototype.merradius = function(B){
	var mr = this.a * (1 - this.e2) * pow((1 - this.e2 * pow(sin (B), 2)), -1.5);
	return mr;
}

datum.prototype.vertradius = function(B){
	var vr = this.a * pow((1 - this.e2 * pow(sin (B), 2)), -0.5);
	return vr;
}

function proj(pr_set){
	this.Name = pr_set.pname;
	this.Ptype = pr_set.ptype;
	this.z = pr_set.pz;
	this.S = pr_set.pscale;
	this.B0 = pr_set.pB0 * PI / 180;
	this.L0 = pr_set.pL0 * PI / 180;
	this.N0 = pr_set.pN0;
	this.E0 = pr_set.pE0;
}

function BLHtoXYZ(e, BLH){
	var XYZ = {};
	var B = toRadians(BLH.B);
	var L = toRadians(BLH.L);
	var H = BLH.H;
	var N = e.vertradius(B);
	XYZ.X = (N + H) * cos(B) * cos(L);
	XYZ.Y = (N + H) * cos(B) * sin(L);
	XYZ.Z = (N + H - e.e2 * N) * sin(B);
	return XYZ;
}

function XYZtoBLH(e, XYZ){
	var BLH = {};
	var N=0;
	var L = atan(XYZ.Y / XYZ.X);
	var Q = sqrt(pow (XYZ.X, 2) + pow(XYZ.Y, 2));
	var tanB = XYZ.Z / Q * (1 / (1 - e.e2));
	var B = atan(tanB);
	var i;
	for (i = 0; i < 3; i++) {
		N = e.vertradius(B);
		var T = XYZ.Z + N * e.e2 * sin(B);
		B = atan(T / Q);
	}
	var H = Q * cos(B) + XYZ.Z * sin(B) - N * (1 - e.e2 * pow(sin(B), 2));
	BLH.B = toDegrees(B);
	BLH.L = toDegrees(L);
	BLH.H = H;
	return BLH;
}

function BLHtoNEH(e, p, BLH){
	var NEH = {};
	var B = toRadians(BLH.B);
	var L = toRadians(BLH.L);
	var H = BLH.H;
	var X = e.merdist(B);
	var X0 = e.merdist(p.B0);
	var M = e.merradius(B);
	var N = e.vertradius(B);
	var n2 = N / M;
	var l = L - p.L0;
	var sinB = sin (B);
	var cosB = cos (B);
	var tanB = tan (B);
	var TN1 = pow(l, 2) / 2 * N * sinB * cosB;
	var TN2 = pow(l, 4) / 24 * N * sinB * pow(cosB, 3) * (4 * pow(n2, 2) + n2 - pow(tanB, 2));
	var TN3 = pow(l, 6) / 720 * N * sinB * pow(cosB, 5) * (8 * pow(n2, 4) * (11 - 24 * pow(tanB, 2)) - 28 * pow(n2, 3) * (1 - 6 * pow(tanB, 2))+ pow(n2, 2) * (1 - 32 * pow(tanB, 2)) -  n2 * (2 * pow(tanB, 2)) + pow(tanB, 4));
	var TN4 = pow(l, 8) / 40320 * N * sinB * pow(cosB, 7) * B * (1385 - 3111 * pow(tanB, 2) + 543 * pow(tanB, 4) - pow(tanB, 6));
	NEH.N = p.N0 + p.S * (X - X0 + TN1 + TN2 + TN3 + TN4);
	var TE1 = pow(l, 2) / 6 * pow(cosB, 2) * (n2 - pow(tanB, 2));
	var TE2 = pow(l, 4) / 120 * pow(cosB, 4) * (4 * pow(n2, 3) * (1 - 6 * pow(tanB, 2)) + pow(n2, 2) * (1 + 8 * pow(tanB, 2)) - n2 * 2 * pow(tanB, 2) + pow(tanB, 4));
	var TE3 = pow(l, 6) / 5040 * pow(cosB, 6) * (61 - 479 * pow(tanB, 2) + 179 * pow(tanB, 4) - pow(tanB, 6));
	NEH.E = p.E0 + p.S * N * l * cosB * (1 + TE1 + TE2 + TE3);
	NEH.H = H;
	return NEH;
}

function NEHtoBLH(e, p, NEH){
	var BLH = {};
	var North = NEH.N;
	var East = NEH.E;
	var X0 = e.merdist(p.B0);
	var N11 = North - p.N0;
	var X1 = X0 + N11 / p.S;
	var n = e.n;
	var q = (X1 * PI) / (180 * e.G);
	var B1 = q + (3 * n / 2 - 27 * pow(n, 3) / 32) * sin (2 * q) + (21 * pow(n, 2) / 16 - 55 * pow(n, 4) / 32) * sin (4 * q) + (151 * pow(n, 3) / 96) * sin (6 * q) + (1097 * pow(n, 4) / 512) * sin (8 * q);
	var M1 = e.merradius(B1);
	var N1 = e.vertradius(B1);
	var n21 = N1 / M1;
	var sinB1 = sin (B1);
	var cosB1 = cos (B1);
	var tanB1 = tan (B1);
	var E1 = East - p.E0;
	var x = E1 / (p.S * N1);
	var TN1 = (tanB1 / (p.S * M1)) * ((E1 * x) / 2);
	var TN2 = (tanB1 / (p.S * M1)) * ((E1 * pow(x, 3)) / 24) * (- 4 * pow(n21, 2) + 9 * n21 * (1 - pow(tanB1, 2)) + 12 * pow(tanB1, 2));
	var TN3 = (tanB1 / (p.S * M1)) * ((E1 * pow(x, 5)) / 720) * (8 * pow(n21, 4) * (11 - 24 * pow(tanB1, 2)) - 12 * pow(n21, 3) * (21 - 71 * pow(tanB1, 2)) + 15 * pow(n21, 2) * (15 - 98 * pow(tanB1, 2) + 15 * pow(tanB1, 4)) + 180 * n21 * (5 * pow(tanB1, 2) - 3 * pow(tanB1, 4)) + 360 * pow(tanB1, 4));
	var TN4 = (tanB1 / (p.S * M1)) * ((E1 * pow(x, 7)) / 40320) * (1385 - 3633 * pow(tanB1, 2) + 4095 * pow(tanB1, 4) + 1575 * pow(tanB1, 6));
	BLH.B = toDegrees(B1 - TN1 + TN2 - TN3 + TN4);
	var TE1 = x * sec(B1);
	var TE2 = pow(x, 3) * sec(B1) / 6 * (n21 + 2 * pow(tanB1, 2));
	var TE3 = pow(x, 5) * sec(B1) / 120 * (-4 * n21 * 3 * (1 - 6 * pow(tanB1, 2)) + pow(n21, 2) * (9 - 68 * pow(tanB1, 2)) + 72 * n21 * pow(tanB1, 2) + 24 * pow(tanB1, 4));
	var TE4 = pow(x, 7) * sec(B1) / 5040 * (61 + 662 * pow(tanB1, 2) + 1320 * pow(tanB1, 4) + 720 * pow(tanB1, 6));
	BLH.L = toDegrees(p.L0 + TE1 - TE2 + TE3 - TE4);
	BLH.H = NEH.H;
	return BLH;
}

 //Position Vector convention of Helmert Transformation Method (Bursa_Wolf Model) (Topcon Sokkia)
 function helmert_pv(e, xyz) {
 	var xyz2 = {};
 	var m = 1 + e.scale/1000000.0;
 	var sinRx = sin(toRadians(e.wx / 3600));
 	var cosRx = cos(toRadians(e.wx / 3600));
 	var sinRy = sin(toRadians(e.wy / 3600));
 	var cosRy = cos(toRadians(e.wy / 3600));
 	var sinRz = sin(toRadians(e.wz / 3600));
 	var cosRz = cos(toRadians(e.wz / 3600));

 	xyz2.X = e.dx + m * (cosRy*cosRz*xyz.X + (-cosRx*sinRz+sinRx*sinRy*cosRz)*xyz.Y + (sinRx*sinRz+cosRx*sinRy*cosRz)*xyz.Z);
 	xyz2.Y = e.dy + m * (cosRy*sinRz*xyz.X + (cosRx*cosRz+sinRx*sinRy*sinRz)*xyz.Y + (-sinRx*cosRz+cosRx*sinRy*sinRz)*xyz.Z);
 	xyz2.Z = e.dz + m * (-sinRy*xyz.X + sinRx*cosRy*xyz.Y + cosRx*cosRy*xyz.Z);  
 	return xyz2
 }

//Coordinate Frame convention of Helmert Transformation Method (Bursa_Wolf Model) (Leica Trimble SP)
 function helmert_cf(e, xyz){
 	var e_cf = e;
 	e_cf.wx = -e.wx;
 	e_cf.wy = -e.wy;
 	e_cf.wz = -e.wz;
 	return helmert_pv(e_cf, xyz);
 }
