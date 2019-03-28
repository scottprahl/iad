\def\title{iad}
\def\iadversion{3-11-1}
\def\iadyear{2019}

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


\def\title{IAD (v \iadversion)}
\def\topofcontents{\null\vfill
  \centerline{\titlefont Inverse Adding-Doubling}
  \vskip 15pt
  \centerline{(Version \iadversion)}
  \vfill}
\def\botofcontents{\vfill
\noindent
Copyright \copyright\ \iadyear\ Scott Prahl
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

%\pageno=\contentspagenumber \advance\pageno by 1

@i "iad_main.w"
@i "iad_type.w"
@i "iad_pub.w"
@i "iad_io.w"
@i "iad_calc.w"
@i "iad_find.w"
@i "iad_util.w"

@**Index.
Here is a cross-reference table for the inverse adding-doubling program.
All sections in which an identifier is
used are listed with that identifier, except that reserved words are
indexed only when they appear in format definitions, and the appearances
of identifiers in section names are not indexed. Underlined entries
correspond to where the identifier was declared. Error messages and
a few other things like ``ASCII code dependencies'' are indexed here too.
