@ECHO OFF
ECHO BASIC VALIDATION TESTS FOR iad.exe PROGRAM
ECHO.

ECHO ********* Basic tests ***********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0
ECHO EXPECT	   0.0000	   0.0000	   0.0000	   0.0000	   1.0000	   0.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 1
ECHO EXPECT	   1.0000	   0.9975	   0.0000	   0.0000	   0.0000	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.1217	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4671	   3.9307	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4346	   3.9791	   0.3116
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.049787
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.6221	   3.9698	  -0.6695
ECHO.

ECHO ********* Specify sample index ************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.0407	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2283	   5.6944	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2274	   5.6964	   0.0354
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4361	   4.5201	  -0.7630
ECHO.

ECHO ********* Specify slide index ************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.0501	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2596	   5.2724	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2564	   5.2748	   0.1021
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4663	   4.4382	  -0.7532
ECHO.

ECHO ********* One slide on top ************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4 -n 1.4 -N 1.5 -G t
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.0501	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -n 1.4 -N 1.5 -G t
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2608	   5.3013	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.4 -N 1.5 -G t
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2577	   5.3037	   0.0991
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5 -G t
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4688	   4.4589	  -0.7536
ECHO.

ECHO ********* One slide on bottom ************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4 -n 1.4 -N 1.5 -G b
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.0487	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -n 1.4 -N 1.5 -G b
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2564	   5.3589	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.4 -N 1.5 -G b
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2536	   5.3618	   0.0898
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5 -G b
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4800	   4.4960	  -0.7760
ECHO.

ECHO ********* Absorbing Slide Tests ***********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.0000000 -t 0.135335 -E 0.5
ECHO EXPECT	   0.0000	   0.0000	   0.1353	   0.1353	   1.0000	   0.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.0249268 -t 0.155858 -E 0.5
ECHO EXPECT	   0.0249	   0.0251	   0.1559	   0.1331	   0.5000	   0.5000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.0520462 -t 0.134587 -E 0.5 -n 1.5 -N 1.5
ECHO EXPECT	   0.0520	   0.0520	   0.1346	   0.1346	   0.5018	   0.4981	   0.0000
ECHO.

ECHO ********* Constrain g ************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4        -g 0.9
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.0101	   0.1000	   0.9000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -g 0.9
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4000	   4.0750	   0.9000
ECHO.

..\iad.exe -V 0 -r 0.4        -g 0.9 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.0040	   0.1000	   0.9000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -g 0.9 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2192	   5.6242	   0.9000
ECHO.

..\iad.exe -V 0 -r 0.4        -g 0.9 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.0048	   0.1000	   0.9000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -g 0.9 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2464	   5.2219	   0.9000
ECHO.

ECHO ********* Constrain a ************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4        -a 0.9
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.1111	   0.9316	   0.0684
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -a 0.9
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4608	   3.9379	   0.0505
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -a 0.9 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.3298	   4.8687	  -0.6404
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -a 0.9 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.3396	   4.7810	  -0.5644
ECHO.

ECHO ********* Constrain b ************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4 -t 0.1 -b 3
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.6221	   3.9698	  -0.6694
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -b 3 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4361	   4.5203	  -0.7631
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -b 3 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4642	   4.4406	  -0.7512
ECHO.

ECHO ********* Constrain mu_s ************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4        -F 30
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   3.6512	  30.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -F 30
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4015	   4.0691	   0.8644
ECHO.

..\iad.exe -V 0 -r 0.4        -F 30 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   1.2216	  30.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -F 30 -n 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2193	   5.6405	   0.8120
ECHO.

..\iad.exe -V 0 -r 0.4        -F 30 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   1.5044	  30.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -F 30 -n 1.4 -N 1.5
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.2467	   5.2304	   0.8257
ECHO.

ECHO ********* Constrain mu_a ************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.3        -A 0.6
ECHO EXPECT	   0.3000	   0.3000	   0.0000	   0.0000	   0.6000	   2.6427	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.3 -t 0.1 -A 0.6
ECHO EXPECT	   0.3000	   0.3000	   0.1000	   0.1000	   0.6000	   3.0846	   0.4187
ECHO.

..\iad.exe -V 0 -r 0.3        -A 0.6 -n 1.5
ECHO EXPECT	   0.3000	   0.3000	   0.0000	   0.0000	   0.6000	   7.2956	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.3 -t 0.1 -A 0.6 -n 1.5
ECHO EXPECT	   0.3000	   0.3000	   0.1000	   0.1000	   0.6000	   3.4619	  -0.6482
ECHO.

..\iad.exe -V 0 -r 0.3        -A 0.6 -n 1.4 -N 1.5
ECHO EXPECT	   0.3000	   0.3000	   0.0000	   0.0000	   0.6000	   5.9837	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.3 -t 0.1 -A 0.6 -n 1.4 -N 1.5
ECHO EXPECT	   0.3000	   0.3000	   0.1000	   0.1000	   0.6000	   3.4104	  -0.5899
ECHO.

ECHO ********* Constrain mu_a and g************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.3        -A 0.6 -g 0.6
ECHO EXPECT	   0.3000	   0.3000	   0.0000	   0.0000	   0.6000	   3.1693	   0.6000
ECHO.

..\iad.exe -V 0 -r 0.3        -A 0.6 -g 0.6 -n 1.5
ECHO EXPECT	   0.3000	   0.3000	   0.0000	   0.0000	   0.6000	   7.7464	   0.6000
ECHO.

..\iad.exe -V 0 -r 0.3        -A 0.6 -g 0.6 -n 1.4 -N 1.5
ECHO EXPECT	   0.3000	   0.3000	   0.0000	   0.0000	   0.6000	   6.4312	   0.6000
ECHO.

ECHO ********* Constrain mu_s and g************
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.3        -F 2.0 -g 0.5
ECHO EXPECT	   0.3000	   0.3000	   0.0000	   0.0000	   0.1939	   1.0000	   0.5000
ECHO.

..\iad.exe -V 0 -r 0.3        -F 2.0 -g 0.5 -n 1.5
ECHO EXPECT	   0.3000	   0.3000	   0.0000	   0.0000	   0.0778	   1.0000	   0.5000
ECHO.

..\iad.exe -V 0 -r 0.3        -F 2.0 -g 0.5 -n 1.4 -N 1.5
ECHO EXPECT	   0.3000	   0.3000	   0.0000	   0.0000	   0.0939	   1.0000	   0.5000
ECHO.

ECHO ********* Basic One Sphere tests ***********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4	   	    -S 1
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.1217	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1	       -S 1
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4671	   3.9307	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002     -S 1
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4346	   3.9791	   0.3116
ECHO.

..\iad.exe -V 0 -r 0.2 -t 0.2 -u 0.0049787 -S 1
ECHO EXPECT	   0.2000	   0.2000	   0.2000	   0.2000	   0.5704	   1.6886	   0.6432
ECHO.

ECHO ******** Basic 10,000 photon tests *********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4	   	    -S 1 -p 10000
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.1217	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1	       -S 1 -p 10000
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4671	   3.9307	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002     -S 1 -p 10000
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4346	   3.9791	   0.3116
ECHO.

..\iad.exe -V 0 -r 0.2 -t 0.2 -u 0.0049787 -S 1 -p 10000
ECHO EXPECT	   0.2000	   0.2000	   0.2000	   0.2000	   0.5704	   1.6886	   0.6432
ECHO.

ECHO ******** Basic timed photon tests *********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4	   	    -S 1 -p -500
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.1217	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1	       -S 1 -p -500
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4671	   3.9307	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002     -S 1 -p -500
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4346	   3.9791	   0.3116
ECHO.

..\iad.exe -V 0 -r 0.2 -t 0.2 -u 0.0049787 -S 1 -p -500
ECHO EXPECT	   0.2000	   0.2000	   0.2000	   0.2000	   0.5704	   1.6886	   0.6432
ECHO.

ECHO ********* More One Sphere tests ***********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4	   	    -S 1 -1 "200 13 13 0"
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.1217	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1	       -S 1 -1 "200 13 13 0"
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4671	   3.9307	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002     -S 1 -1 "200 13 13 0"
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4346	   3.9791	   0.3116
ECHO.

..\iad.exe -V 0 -r 0.2 -t 0.2 -u 0.0049787 -S 1 -1 "200 13 13 0"
ECHO EXPECT	   0.2000	   0.2000	   0.2000	   0.2000	   0.5704	   1.6886	   0.6432
ECHO.

ECHO ********* Basic Two Sphere tests ***********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4	   	    -S 2
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.1217	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1	       -S 2
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4671	   3.9307	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002     -S 2
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4346	   3.9791	   0.3116
ECHO.

..\iad.exe -V 0 -r 0.2 -t 0.1 -u 0.0049787 -S 2
ECHO EXPECT	   0.2000	   0.2000	   0.1000	   0.1000	   0.8339	   2.2615	   0.4939
ECHO.

ECHO ********* More Two Sphere tests ***********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4	   	    -S 2 -1 "200 13 13 0" -2 "200 13 13 0"
ECHO EXPECT	   0.4000	   0.4000	   0.0000	   0.0000	   0.1217	   1.0000	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1	       -S 2 -1 "200 13 13 0" -2 "200 13 13 0"
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4671	   3.9307	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.1 -u 0.002     -S 2 -1 "200 13 13 0" -2 "200 13 13 0"
ECHO EXPECT	   0.4000	   0.4000	   0.1000	   0.1000	   0.4346	   3.9791	   0.3116
ECHO.

..\iad.exe -V 0 -r 0.2 -t 0.1 -u 0.0049787 -S 2 -1 "200 13 13 0" -2 "200 13 13 0"
ECHO EXPECT	   0.2000	   0.2000	   0.1000	   0.1000	   0.8339	   2.2615	   0.4939
ECHO.

ECHO ********* Oblique tests ***********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g
f
..\iad.exe -V 0 -i 60 -r 0.00000 -t 0.13691
ECHO EXPECT	   0.0000	   0.0000	   0.1369	   0.1369	   0.9941	   0.0000	   0.0000
ECHO.

..\iad.exe -V 0 -i 60 -r 0.14932 -t 0.23181
ECHO EXPECT	   0.1493	   0.1493	   0.2318	   0.2318	   0.4980	   0.4966	   0.0000
ECHO.

..\iad.exe -V 0 -i 60 -r 0.61996 -t 0.30605
ECHO EXPECT	   0.6200	   0.6199	   0.3060	   0.3060	   0.0382	   2.0176	   0.0000
ECHO.

ECHO ********* Different Tstd and Rstd ***********
ECHO.

ECHO 	     Meas R	   Calc R	   Meas T	   Calc T	     mu_a	    mu_s"	        g

..\iad.exe -V 0 -r 0.4 -t 0.005  -d 1 -M 0  -S 1 -1 "200 13 13 0" -T 0.5
ECHO EXPECT	   0.4000	   0.4000	   0.0050	   0.0050	   1.0739	   8.8255	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.4 -t 0.005  -d 1 -M 0  -S 1 -1 "200 13 13 0"
ECHO EXPECT	   0.4000	   0.4000	   0.0050	   0.0050	   1.0739	   8.8255	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.2 -t 0.01   -d 1 -M 0  -S 1 -1 "200 13 13 0" -R 0.5
ECHO EXPECT	   0.2000	   0.2000	   0.0100	   0.0100	   1.8796	   4.1017	   0.0000
ECHO.

..\iad.exe -V 0 -r 0.2 -t 0.01   -d 1 -M 0  -S 1 -1 "200 13 13 0"
ECHO EXPECT	   0.2000	   0.2000	   0.0100	   0.0100	   1.8796	   4.1017	   0.0000
ECHO.
PAUSE
