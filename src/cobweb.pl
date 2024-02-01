#!/usr/bin/perl -pi.bak -w

BEGIN{undef $/;}
s!/\*.*\n!!mg;
s!#line.*\n!!mg;
s!\n\n+!\n\n!sg;
s!\n+;!;!sg;
s!^\n+!!s;
