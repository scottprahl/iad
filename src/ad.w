\def\title{ad}
\def\adversion{3-11-1}
\def\adyear{2019}
% This program by Scott Prahl
% It is distributed WITHOUT ANY WARRANTY, express or implied.

% Copyright (C) 2019 Scott Prahl

% Permission is granted to make and distribute verbatim copies of this
% document provided that the copyright notice and this permission notice
% are preserved on all copies.

% Permission is granted to copy and distribute modified versions of this
% document under the conditions for verbatim copying, provided that the
% entire resulting derived work is given a different name and distributed
% under the terms of a permission notice identical to this one.

\def\title{Adding-Doubling (Version \adversion)}
\def\topofcontents{\null\vfill
  \centerline{\titlefont The Adding-Doubling Program}
  \vskip 15pt
  \centerline{(Version \adversion)}
  \vfill}
\def\botofcontents{\vfill
\noindent
Copyright \copyright\ \adyear\ Scott Prahl
\bigskip\noindent
Permission is granted to make and distribute verbatim copies of this
document provided that the copyright notice and this permission notice
are preserved on all copies.

\smallskip\noindent
Permission is granted to copy and distribute modified versions of this
document under the conditions for verbatim copying, provided that the
entire resulting derived work is given a different name and distributed
under the terms of a permission notice identical to this one.
}

@i "ad_globl.w"
@i "ad_prime.w"
@i "ad_layers.w"
@i "ad_cone.w"
@i "ad_start.w"
@i "ad_doubl.w"
@i "ad_bound.w"
@i "ad_frsnl.w"
@i "ad_matrx.w"
@i "ad_radau.w"
@i "ad_phase.w"
@i "ad_main.w"

@**Index.
Here is a cross-reference table for the adding-doubling program.
All sections in which an identifier is
used are listed with that identifier, except that reserved words are
indexed only when they appear in format definitions, and the appearances
of identifiers in section names are not indexed. Underlined entries
correspond to where the identifier was declared. Error messages and
a few other things like ``ASCII code dependencies'' are indexed here too.
