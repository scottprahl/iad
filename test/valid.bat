@ECHO OFF
ECHO BASIC VALIDATION TESTS FOR iad.exe PROGRAM
ECHO.

ECHO ********* Basic tests ***********
..\iad.exe -V0 -r 0
ECHO EXPECT  0.0000  0.0000  0.0000  0.0000  0.0000  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 1
ECHO EXPECT  1.0000  0.9975  0.0000  0.0000  0.0000  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.1217  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.4671  3.9307  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.01
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.4609  3.9375  0.0499
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.049787
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.6221  3.9698  -0.6695
ECHO.

ECHO.
ECHO ********* Specify sample index ************
..\iad.exe -V0 -r 0.4 -n 1.5
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0407  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -n 1.5
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2283  5.6944  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.01 -n 1.5
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2439  5.5666  -0.3007
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.045884 -n 1.5
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.4361  4.5201  -0.7630
ECHO.

ECHO.
ECHO ********* Specify slide index ************
..\iad.exe -V0 -r 0.4 -n 1.4 -N 1.5
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0501  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -n 1.4 -N 1.5
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2596  5.2724  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.01 -n 1.4 -N 1.5
ECHO EXPECT  0.4000  0.4221  0.1000  0.1000  0.8929  7.2568  -1.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.6197  4.4494  -0.8710
ECHO.

ECHO.
ECHO ********* Constrain g ************
..\iad.exe -V0 -r 0.4        -g 0.9
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0101  0.1000  0.9000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -g 0.9
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.4000  4.0750  0.9000
ECHO.

..\iad.exe -V0 -r 0.4        -g 0.9 -n 1.5
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0040  0.1000  0.9000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -g 0.9 -n 1.5
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2192  5.6242  0.9000
ECHO.

..\iad.exe -V0 -r 0.4        -g 0.9 -n 1.4 -N 1.5
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0048  0.1000  0.9000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -g 0.9 -n 1.4 -N 1.5
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2464  5.2219  0.9000
ECHO.

ECHO.
ECHO ********* Constrain a ************
..\iad.exe -V0 -r 0.4        -a 0.9
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.1111  0.9316  0.0684
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -a 0.9
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  4.6080  0.0000  0.0504
ECHO.

ECHO.
ECHO ********* Constrain b ************
..\iad.exe -V0 -r 0.4        -b 3
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.1217  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -b 3
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.6221  3.9698  -0.6695
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -b 3 -n 1.5
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.4361  4.5201  -0.7630
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -b 3 -n 1.4 -N 1.5
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.6223  4.4527  -0.8727
ECHO.

ECHO.
ECHO ********* Basic One Sphere tests ***********
..\iad.exe -V0 -r 0.4                     -S 1
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0482  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1              -S 1
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2683  6.8248  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.01      -S 1
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.3131  6.9026  -0.6082
ECHO.

..\iad.exe -V0 -r 0.2 -t 0.2 -u 0.049787  -S 1
ECHO EXPECT  0.2000  0.2000  0.2000  0.2000  0.4665  2.5940  -0.0239
ECHO.

ECHO.
ECHO ******** Basic 10,000 photon tests *********
..\iad.exe -V0 -r 0.4                     -S 1 -p 10000
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0482  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1              -S 1 -p 10000
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2704  6.7994  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.01      -S 1 -p 10000
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.3137  6.8921  -0.6060
ECHO.

..\iad.exe -V0 -r 0.2 -t 0.2 -u 0.049787  -S 1 -p 10000
ECHO EXPECT  0.2000  0.2000  0.2000  0.2000  0.4659  2.5955  -0.0242
ECHO.

ECHO.
ECHO ******** Basic timed photon tests *********
..\iad.exe -V0 -r 0.4                     -S 1 -p -100
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0482  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1              -S 1 -p -100
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2698  6.7850  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.01      -S 1 -p -100
ECHO EXPECT  0.4000  0.3999  0.1000  0.1000  0.3111  6.9382  -0.6158
ECHO.

..\iad.exe -V0 -r 0.2 -t 0.2 -u 0.049787  -S 1 -p -100
ECHO EXPECT  0.2000  0.2000  0.2000  0.2000  0.4670  2.5914  -0.0231
ECHO.

ECHO.
ECHO ********* More One Sphere tests ***********
..\iad.exe -V0 -r 0.4                     -S 1 -1 '200 13 13 0'
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0482  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1              -S 1 -1 '200 13 13 0'
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2683  6.8248  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.01      -S 1 -1 '200 13 13 0'
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.3131  6.9026  -0.6082
ECHO.

..\iad.exe -V0 -r 0.2 -t 0.2 -u 0.049787  -S 1 -1 '200 13 13 0'
ECHO EXPECT  0.2000  0.2000  0.2000  0.2000  0.4665  2.5940  -0.0239
ECHO.

ECHO.
ECHO ********* Basic Two Sphere tests ***********
..\iad.exe -V0 -r 0.4                     -S 2
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0482  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1              -S 2
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2657  6.6591  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.01      -S 2
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.3004  6.6607  -0.5473
ECHO.

..\iad.exe -V0 -r 0.2 -t 0.1 -u 0.049787  -S 2
ECHO EXPECT  0.2000  0.2000  0.1000  0.1000  0.7048  3.2091  -0.3982
ECHO.

ECHO.
ECHO ********* More Two Sphere tests ***********
..\iad.exe -V0 -r 0.4                     -S 2 -1 '200 13 13 0' -2 '200 13 13 0'
ECHO EXPECT  0.4000  0.4000  0.0000  0.0000  0.0482  1.0000  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1              -S 2 -1 '200 13 13 0' -2 '200 13 13 0'
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.2657  6.6591  0.0000
ECHO.

..\iad.exe -V0 -r 0.4 -t 0.1 -u 0.01      -S 2 -1 '200 13 13 0' -2 '200 13 13 0'
ECHO EXPECT  0.4000  0.4000  0.1000  0.1000  0.3004  6.6607  -0.5473
ECHO.

..\iad.exe -V0 -r 0.2 -t 0.1 -u 0.049787  -S 2 -1 '200 13 13 0' -2 '200 13 13 0'
ECHO EXPECT  0.2000  0.2000  0.1000  0.1000  0.7048  3.2091  -0.3982
ECHO.
PAUSE