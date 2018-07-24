%%
%% This is the file `afmc_art.cls'
%% generated manually by P.V. Lanspeary for
%% Australasian Fluid Mechanics Conference papers
%% (except where it is copied from 13AFMC stencil.tex & emac98.sty :-)
%% Jason Roberts : redefine ``thebibliography'' for  to fix bug with Reference label.
%% TL made some tiny adjustments to PVL's original file for 21AFMC

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{afmc_art_2018}[2018/03/01 AFMC article class]

%% The AFMC template is a very narrow specification:
\DeclareOption{onecolumn}{\OptionNotUsed}
\DeclareOption{11pt}{\OptionNotUsed}
\DeclareOption{12pt}{\OptionNotUsed}
\DeclareOption{letterpaper}{\OptionNotUsed}
\DeclareOption{legalpaper}{\OptionNotUsed}
\DeclareOption{executivepaper}{\OptionNotUsed}
\DeclareOption{a5paper}{\OptionNotUsed}
\DeclareOption{b5paper}{\OptionNotUsed}
\DeclareOption{titlepage}{\OptionNotUsed}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions

\LoadClass[a4paper,10pt,twocolumn,notitlepage,twoside]{article}
\usepackage{pslatex}
%%
%% Page Layout
%%
\renewcommand{\columnsep}{10mm}
\setlength{\textheight}{248mm}
\setlength{\textwidth}{170mm}
\setlength{\oddsidemargin}{-5.4mm}    % 30 mm left margin on 1st page
\setlength{\evensidemargin}{-5.4mm}   % 20 mm left margin on 2nd page
\setlength{\topmargin}{5mm}
\setlength{\headsep}{-10mm}

%%
%% Font sizes
%% use of the cmr/cms bitmap fonts is discontinued because we
%% are using the pslatex fonts -- for the sake of nice pdf
%%
\renewcommand{\normalsize}{\fontsize{9}{9}\selectfont}
\newcommand{\authorsize}{\fontsize{10}{10}\selectfont}
\newcommand{\titlesize}{\fontsize{12}{12}\selectfont}

%%
%% Sections etc
%% --- Just mimic what was in the original template, nothing fancy
%%     If we try to simply modify the definitions in article.cls,
%%     we find that we cannot reduce the space between the section
%%     title and the first line to zero ... hence the definitions
%%     below - which could be done equally well on a typewriter !
%%     Of course, without numbering, one cannot label and use \ref.
%%                                 }
\renewcommand{\section}[1]{\vspace{1.0ex}%
                           {\sffamily \bfseries \normalsize {#1}}%
                           \par \vspace{-0.5ex}}
\renewcommand{\subsection}[1]{\vspace{1.0ex}%
                              {\sffamily \normalsize \underline{#1}}
                              \par \vspace{-0.5ex}}
\renewcommand{\subsubsection}[1]{\typeout{Please do not use subsubsection!}}
\renewcommand{\paragraph}[1]{\typeout{Please do not use paragraph!}\par}
\renewcommand{\subparagraph}[1]{\typeout{Please do not use subparagraph!}\par}

%% Abstract:
%% use of unnumbered sections and two-column page layout
%% makes the abstract just another section
\renewenvironment{abstract}{\section{Abstract}}{}
\setcounter{secnumdepth}{0}

%%
%% Page headings
%%  The first page contains the heading
%%    "21st Australasian Fluid Mechanics Conference
%%     Adelaide, Australia
%%     10-13 December 2018"
%%  So we define a new page style called afmchead and
%%  invoke it in the \maketitle.  The incantation is
%%  a modification of that found in Kopka & Daly, Appendix C.3.2
%%  The other pages contain no heading.

\pagestyle{empty}

\newcommand{\ps@afmchead}{%
  \renewcommand{\@oddhead}{\setlength{\unitlength}{1pt}%
    \renewcommand{\baselinestretch}{1.00}%
    \begin{picture}(-0,0)(0,0)
    \put(0,40){\parbox{100mm}{\rmfamily \mdseries \normalsize
         21st Australasian Fluid Mechanics Conference\\
         Adelaide, Australia\\
         10-13 December 2018}}
    \end{picture}}
  \renewcommand{\@evenhead}{}
  \renewcommand{\@oddfoot}{}
  \renewcommand{\@evenfoot}{}
}


\raggedbottom
%
% define paragraph format
\setlength{\parindent}{0.0em}
\setlength{\parskip}{1.5ex}
\setlength{\textfloatsep}{16pt plus 2pt minus 4pt}
%
%  expand the line spacing a bit. Is this enough ?
%  In future versions, it might be better use \fontsize for this?
\renewcommand{\baselinestretch}{1.15}
%


%
% set-up definitions for title (pinched from emac98.sty)
%
\def\affiliation#1{\gdef\@affiliation{#1}}
\gdef\@affiliation{}
\renewcommand{\maketitle}{\par
  \begingroup
    \renewcommand{\thefootnote}{\fnsymbol{footnote}}%
     \def\@makefnmark{\hbox to\z@{$\m@th^{\@thefnmark}$\hss}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
            \hbox to1.8em{\hss$\m@th^{\@thefnmark}$}##1}%
    \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{afmchead}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \let\thanks\relax
\let\maketitle\relax \let\@maketitle\relax \gdef\@affiliation{}
 \gdef\@thanks{}\gdef\@author{}\gdef\@title{}\let\thanks\relax}
\def\@maketitle{%
  \newpage
 \null
 \begin{center}%
%  \vskip -2.5em%
   {\sffamily \bfseries \titlesize \@title \par}%
   \vskip 1.5em%
   {\sffamily \bfseries \authorsize \@author \par}%
   \vskip 1.5em%
   {\sffamily \mdseries \authorsize \@affiliation \par}%
 \end{center}%
 \par
 \vskip 1.0em}

\def\refname{\sffamily \bfseries \normalsize References}
%
% Addition by Jason Roberts <J.L.Roberts@utas.edu.au>
% to fix bug with heading of the References section.
%
\renewenvironment{thebibliography}[1]
     {\vspace{1.0ex}%
      \section{\refname
        \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
      \par \vspace{-0.5ex}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}

%
% ignore requests for tableofcontents, listoffigures, listoftables.
% Discourage index.
%
\renewcommand{\tableofcontents}
             {\typeout{Request for tableofcontents is ignored.}}
\renewcommand{\listoffigures}
             {\typeout{Request for listoffigures is ignored.}}
\renewcommand{\listoftables}
             {\typeout{Request for listoftables is ignored.}}
\renewenvironment{theindex}
               {\typeout{Please do not put an index in your AFMC paper!}
                \if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\section{\indexname}]%
                \@mkboth{\MakeUppercase\indexname}%
                        {\MakeUppercase\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
