/*2:*/
#line 32 "ad_doubl.w"

/*6:*/
#line 184 "ad_doubl.w"

void Add(int n,
double**R01,double**R10,double**T01,double**T10,
double**R12,double**R21,double**T12,double**T21,
double**R02,double**R20,double**T02,double**T20)

/*:6*/
#line 33 "ad_doubl.w"
;
/*8:*/
#line 208 "ad_doubl.w"

void Add_With_Sources(int n,
double**R01,double**R10,double**T01,double**T10,double**J01,double**J10,
double**R12,double**R21,double**T12,double**T21,double**J12,double**J21,
double**R02,double**R20,double**T02,double**T20,double**J02,double**J20)

/*:8*/
#line 34 "ad_doubl.w"
;
/*10:*/
#line 232 "ad_doubl.w"

void Add_Homogeneous(int n,
double**R01,double**T01,
double**R12,double**T12,
double**R02,double**T02)

/*:10*/
#line 35 "ad_doubl.w"
;
/*12:*/
#line 254 "ad_doubl.w"

void Double_Once(int n,double**R,double**T)

/*:12*/
#line 36 "ad_doubl.w"
;
/*14:*/
#line 273 "ad_doubl.w"

void Double_Until(int n,double**r,double**t,double start,double end)

/*:14*/
#line 37 "ad_doubl.w"
;
/*16:*/
#line 305 "ad_doubl.w"

void Double_Until_Infinite(int n,double**r,double**t)

/*:16*/
#line 38 "ad_doubl.w"
;
/*18:*/
#line 352 "ad_doubl.w"

void Between(int n,
double**R01,double**R10,double**T01,double**T10,
double**R12,double**R21,double**T12,double**T21,
double**Lup,double**Ldown)

/*:18*/
#line 39 "ad_doubl.w"
;

/*:2*/
